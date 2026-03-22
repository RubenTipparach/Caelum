using System.IO;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Media3D;
using System.Windows.Threading;
using CharEditor.Generators;
using CharEditor.Models;
using CharEditor.Services;
using Microsoft.Win32;

namespace CharEditor;

public partial class MainWindow : Window
{
    private readonly AgentCharacter _character = new();
    private readonly ChatService _chat = new();
    private bool _suppressEvents;
    private bool _initialized;

    // Animation state
    private AnimState _animState = AnimState.Idle;
    private HexBodyGenerator.BodyParts? _parts;
    private readonly DispatcherTimer _animTimer;
    private DateTime _animStart = DateTime.Now;
    private bool _isThinking;

    // Animated model nodes
    private GeometryModel3D? _mdlTorso, _mdlHead, _mdlEyes;
    private GeometryModel3D? _mdlLeftArm, _mdlRightArm, _mdlLeftLeg, _mdlRightLeg;

    // Cached transforms (reused every frame to avoid GC pressure)
    private readonly AxisAngleRotation3D _bodyRot = new(new Vector3D(1, 0, 0), 0);
    private readonly TranslateTransform3D _bodyOff = new();
    private readonly AxisAngleRotation3D _headTilt = new(new Vector3D(1, 0, 0), 0);
    private readonly AxisAngleRotation3D _headTurn = new(new Vector3D(0, 1, 0), 0);
    private readonly AxisAngleRotation3D _lArmPitch = new(new Vector3D(1, 0, 0), 0);
    private readonly AxisAngleRotation3D _lArmRaise = new(new Vector3D(0, 0, 1), 0);
    private readonly AxisAngleRotation3D _rArmPitch = new(new Vector3D(1, 0, 0), 0);
    private readonly AxisAngleRotation3D _rArmRaise = new(new Vector3D(0, 0, 1), 0);
    private readonly AxisAngleRotation3D _lLegPitch = new(new Vector3D(1, 0, 0), 0);
    private readonly AxisAngleRotation3D _rLegPitch = new(new Vector3D(1, 0, 0), 0);

    public MainWindow()
    {
        _suppressEvents = true;
        InitializeComponent();

        // Populate combos with plain strings (avoids WPF ComboBoxItem styling issues)
        foreach (var t in DirectiveTemplates.All)
            TemplateCombo.Items.Add(t.Name);
        TemplateCombo.SelectedIndex = 0;

        ClassCombo.Items.Add("Worker");
        ClassCombo.Items.Add("Wanderer");
        ClassCombo.Items.Add("Villager");
        ClassCombo.SelectedIndex = 0;

        _character.Seed = (uint)Random.Shared.Next();
        _suppressEvents = false;
        _initialized = true;
        UpdateValueLabels();
        RebuildMesh();

        // Animation timer — 30fps
        _animTimer = new DispatcherTimer { Interval = TimeSpan.FromMilliseconds(33) };
        _animTimer.Tick += OnAnimTick;
        _animTimer.Start();
        _animStart = DateTime.Now;

        // Chat: check if llama-server is available
        _ = CheckChatServerAsync();
    }

    private async Task CheckChatServerAsync()
    {
        ChatStatus.Text = "Starting llama-server...";
        ChatStatus.Foreground = new SolidColorBrush(Color.FromRgb(0xff, 0xff, 0x66));

        _chat.OnStatusUpdate = msg => Dispatcher.Invoke(() =>
        {
            ChatStatus.Text = msg;
        });

        var result = await _chat.TryAutoLaunchAsync();
        var ok = _chat.IsAvailable;

        ChatStatus.Text = ok ? "llama-server connected" : result;
        ChatStatus.Foreground = new SolidColorBrush(ok
            ? Color.FromRgb(0x66, 0xff, 0x66)
            : Color.FromRgb(0xff, 0x66, 0x66));
    }

    protected override void OnClosed(EventArgs e)
    {
        _animTimer.Stop();
        base.OnClosed(e);
    }

    // ---- Event handlers ----

    private void OnAppearanceChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
    {
        if (_suppressEvents) return;
        SyncFromUi();
        UpdateValueLabels();
        RebuildMesh();
    }

    private void OnPersonalityChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
    {
        if (_suppressEvents) return;
        SyncFromUi();
    }

    private void OnAttributeChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
    {
        if (_suppressEvents) return;
        SyncFromUi();
        UpdateValueLabels();
    }

    private void OnAnyChanged(object sender, object e)
    {
        if (_suppressEvents) return;
        SyncFromUi();
    }

    private void OnTemplateChanged(object sender, SelectionChangedEventArgs e)
    {
        if (_suppressEvents || TemplateCombo.SelectedIndex < 0) return;
        var tmpl = DirectiveTemplates.All[TemplateCombo.SelectedIndex];

        _suppressEvents = true;
        DirectiveBox.Text = tmpl.Directive;
        CreativitySlider.Value = tmpl.Creativity;
        ChattinessSlider.Value = tmpl.Chattiness;
        DiligenceSlider.Value = tmpl.Diligence;
        _suppressEvents = false;

        SyncFromUi();
    }

    private void OnRandomize(object sender, RoutedEventArgs e)
    {
        _character.Randomize(Random.Shared);
        SyncToUi();
        RebuildMesh();
        StatusText.Text = "Randomized!";
    }

    private static string GetAgentsDir()
    {
        // Walk up from the exe to find the project root (has CMakeLists.txt)
        var dir = AppDomain.CurrentDomain.BaseDirectory;
        for (int i = 0; i < 8; i++)
        {
            var parent = Path.GetDirectoryName(dir);
            if (parent == null) break;
            dir = parent;
            if (File.Exists(Path.Combine(dir, "CMakeLists.txt")))
            {
                var agents = Path.Combine(dir, "cache", "agents");
                Directory.CreateDirectory(agents);
                return agents;
            }
        }
        // Fallback: next to the exe
        var fallback = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "agents");
        Directory.CreateDirectory(fallback);
        return fallback;
    }

    private void OnExport(object sender, RoutedEventArgs e)
    {
        SyncFromUi();
        var data = AgentExportData.FromCharacter(_character);
        var json = data.ToJson();

        var agentsDir = GetAgentsDir();
        var dlg = new SaveFileDialog
        {
            Filter = "JSON files (*.json)|*.json",
            FileName = $"{_character.Name.ToLowerInvariant().Replace(' ', '_')}_{_character.Id[..8]}.json",
            Title = "Save Agent Character",
            InitialDirectory = agentsDir
        };
        if (dlg.ShowDialog() == true)
        {
            File.WriteAllText(dlg.FileName, json);
            StatusText.Text = $"Saved to {dlg.FileName}";
        }
    }

    private void OnImport(object sender, RoutedEventArgs e)
    {
        var agentsDir = GetAgentsDir();
        var dlg = new OpenFileDialog
        {
            Filter = "JSON files (*.json)|*.json",
            Title = "Load Agent Character",
            InitialDirectory = agentsDir
        };
        if (dlg.ShowDialog() == true)
        {
            try
            {
                var json = File.ReadAllText(dlg.FileName);
                var data = AgentExportData.FromJson(json);
                if (data == null)
                {
                    StatusText.Text = "Failed to parse JSON.";
                    return;
                }
                var imported = AgentExportData.ToCharacter(data);

                // Copy into our character
                _character.Name = imported.Name;
                _character.Class = imported.Class;
                _character.Seed = imported.Seed;
                _character.Directive = imported.Directive;

                _character.Appearance.Height = imported.Appearance.Height;
                _character.Appearance.Stockiness = imported.Appearance.Stockiness;
                _character.Appearance.HeadScale = imported.Appearance.HeadScale;
                _character.Appearance.ArmLength = imported.Appearance.ArmLength;
                _character.Appearance.TorsoWidth = imported.Appearance.TorsoWidth;
                _character.Appearance.LegLength = imported.Appearance.LegLength;
                _character.Appearance.PrimaryR = imported.Appearance.PrimaryR;
                _character.Appearance.PrimaryG = imported.Appearance.PrimaryG;
                _character.Appearance.PrimaryB = imported.Appearance.PrimaryB;
                _character.Appearance.AccentR = imported.Appearance.AccentR;
                _character.Appearance.AccentG = imported.Appearance.AccentG;
                _character.Appearance.AccentB = imported.Appearance.AccentB;

                _character.Personality.Creativity = imported.Personality.Creativity;
                _character.Personality.Chattiness = imported.Personality.Chattiness;
                _character.Personality.Diligence = imported.Personality.Diligence;

                _character.Attributes.Strength = imported.Attributes.Strength;
                _character.Attributes.Perception = imported.Attributes.Perception;
                _character.Attributes.Endurance = imported.Attributes.Endurance;

                SyncToUi();
                RebuildMesh();
                StatusText.Text = $"Imported {Path.GetFileName(dlg.FileName)}";
            }
            catch (Exception ex)
            {
                StatusText.Text = $"Import error: {ex.Message}";
            }
        }
    }

    // ---- Sync UI <-> Model ----

    private void SyncFromUi()
    {
        if (_suppressEvents || !_initialized) return;

        _character.Name = NameBox.Text;
        _character.Class = (AgentClass)ClassCombo.SelectedIndex;
        _character.Directive = DirectiveBox.Text;

        var a = _character.Appearance;
        a.Height = HeightSlider.Value;
        a.Stockiness = StockySlider.Value;
        a.HeadScale = HeadSlider.Value;
        a.ArmLength = ArmSlider.Value;
        a.TorsoWidth = TorsoSlider.Value;
        a.LegLength = LegSlider.Value;
        a.PrimaryR = PrimR.Value;
        a.PrimaryG = PrimG.Value;
        a.PrimaryB = PrimB.Value;
        a.AccentR = AccR.Value;
        a.AccentG = AccG.Value;
        a.AccentB = AccB.Value;

        _character.Personality.Creativity = CreativitySlider.Value;
        _character.Personality.Chattiness = ChattinessSlider.Value;
        _character.Personality.Diligence = DiligenceSlider.Value;

        _character.Attributes.Strength = (int)StrSlider.Value;
        _character.Attributes.Perception = (int)PerSlider.Value;
        _character.Attributes.Endurance = (int)EndSlider.Value;

        // Update color previews
        PrimPreview.Fill = new SolidColorBrush(
            HexBodyGenerator.ToColor(a.PrimaryR, a.PrimaryG, a.PrimaryB));
        AccPreview.Fill = new SolidColorBrush(
            HexBodyGenerator.ToColor(a.AccentR, a.AccentG, a.AccentB));
    }

    private void SyncToUi()
    {
        _suppressEvents = true;

        NameBox.Text = _character.Name;
        ClassCombo.SelectedIndex = (int)_character.Class;
        DirectiveBox.Text = _character.Directive;

        var a = _character.Appearance;
        HeightSlider.Value = a.Height;
        StockySlider.Value = a.Stockiness;
        HeadSlider.Value = a.HeadScale;
        ArmSlider.Value = a.ArmLength;
        TorsoSlider.Value = a.TorsoWidth;
        LegSlider.Value = a.LegLength;
        PrimR.Value = a.PrimaryR;
        PrimG.Value = a.PrimaryG;
        PrimB.Value = a.PrimaryB;
        AccR.Value = a.AccentR;
        AccG.Value = a.AccentG;
        AccB.Value = a.AccentB;

        CreativitySlider.Value = _character.Personality.Creativity;
        ChattinessSlider.Value = _character.Personality.Chattiness;
        DiligenceSlider.Value = _character.Personality.Diligence;

        StrSlider.Value = _character.Attributes.Strength;
        PerSlider.Value = _character.Attributes.Perception;
        EndSlider.Value = _character.Attributes.Endurance;

        PrimPreview.Fill = new SolidColorBrush(
            HexBodyGenerator.ToColor(a.PrimaryR, a.PrimaryG, a.PrimaryB));
        AccPreview.Fill = new SolidColorBrush(
            HexBodyGenerator.ToColor(a.AccentR, a.AccentG, a.AccentB));

        UpdateValueLabels();
        _suppressEvents = false;
    }

    private void UpdateValueLabels()
    {
        HeightVal.Text = $"{HeightSlider.Value:F1}";
        StockyVal.Text = $"{StockySlider.Value:F1}";
        HeadVal.Text = $"{HeadSlider.Value:F1}";
        ArmVal.Text = $"{ArmSlider.Value:F1}";
        TorsoVal.Text = $"{TorsoSlider.Value:F1}";
        LegVal.Text = $"{LegSlider.Value:F1}";
        StrVal.Text = $"{(int)StrSlider.Value}";
        PerVal.Text = $"{(int)PerSlider.Value}";
        EndVal.Text = $"{(int)EndSlider.Value}";
    }

    // ---- 3D Mesh + Animation ----

    private void RebuildMesh()
    {
        var a = _character.Appearance;
        _parts = HexBodyGenerator.GenerateParts(a);

        var primaryColor = HexBodyGenerator.ToColor(a.PrimaryR, a.PrimaryG, a.PrimaryB);
        var accentColor = HexBodyGenerator.ToColor(a.AccentR, a.AccentG, a.AccentB);

        var primaryMat = new DiffuseMaterial(new SolidColorBrush(primaryColor));
        var accentMat = new DiffuseMaterial(new SolidColorBrush(accentColor));
        var eyeMat = new DiffuseMaterial(new SolidColorBrush(Colors.White));

        _mdlTorso = new GeometryModel3D(_parts.Torso, primaryMat);
        _mdlHead = new GeometryModel3D(_parts.HeadAccent, accentMat);
        _mdlEyes = new GeometryModel3D(_parts.Eyes, eyeMat);
        _mdlLeftArm = new GeometryModel3D(_parts.LeftArm, accentMat);
        _mdlRightArm = new GeometryModel3D(_parts.RightArm, accentMat);
        _mdlLeftLeg = new GeometryModel3D(_parts.LeftLeg, primaryMat);
        _mdlRightLeg = new GeometryModel3D(_parts.RightLeg, primaryMat);

        var group = new Model3DGroup();
        group.Children.Add(_mdlTorso);
        group.Children.Add(_mdlHead);
        group.Children.Add(_mdlEyes);
        group.Children.Add(_mdlLeftArm);
        group.Children.Add(_mdlRightArm);
        group.Children.Add(_mdlLeftLeg);
        group.Children.Add(_mdlRightLeg);

        CharacterModel.Content = group;
        _transformsBuilt = false;
    }

    private void OnAnimSelect(object sender, RoutedEventArgs e)
    {
        if (sender is Button btn && btn.Tag is string tag)
        {
            if (Enum.TryParse<AnimState>(tag, out var state))
            {
                _animState = state;
                _animStart = DateTime.Now;
            }
        }
    }

    private bool _transformsBuilt;

    private void BuildTransforms()
    {
        if (_parts == null || _transformsBuilt) return;

        var feetPivot = new Point3D(0, 0, 0);
        var headPivot = new Point3D(0, _parts.HeadPivotY, 0);
        var lArmPivot = new Point3D(-_parts.ArmSpacing, _parts.ArmPivotY, 0);
        var rArmPivot = new Point3D(_parts.ArmSpacing, _parts.ArmPivotY, 0);
        var lLegPivot = new Point3D(-_parts.LegSpacing, _parts.LegPivotY, 0);
        var rLegPivot = new Point3D(_parts.LegSpacing, _parts.LegPivotY, 0);

        var bodyRot = new RotateTransform3D(_bodyRot, feetPivot);

        var torsoTf = new Transform3DGroup();
        torsoTf.Children.Add(bodyRot);
        torsoTf.Children.Add(_bodyOff);
        _mdlTorso!.Transform = torsoTf;

        var headTf = new Transform3DGroup();
        headTf.Children.Add(new RotateTransform3D(_headTilt, headPivot));
        headTf.Children.Add(new RotateTransform3D(_headTurn, headPivot));
        headTf.Children.Add(bodyRot);
        headTf.Children.Add(_bodyOff);
        _mdlHead!.Transform = headTf;
        _mdlEyes!.Transform = headTf;

        var lArmTf = new Transform3DGroup();
        lArmTf.Children.Add(new RotateTransform3D(_lArmPitch, lArmPivot));
        lArmTf.Children.Add(new RotateTransform3D(_lArmRaise, lArmPivot));
        lArmTf.Children.Add(bodyRot);
        lArmTf.Children.Add(_bodyOff);
        _mdlLeftArm!.Transform = lArmTf;

        var rArmTf = new Transform3DGroup();
        rArmTf.Children.Add(new RotateTransform3D(_rArmPitch, rArmPivot));
        rArmTf.Children.Add(new RotateTransform3D(_rArmRaise, rArmPivot));
        rArmTf.Children.Add(bodyRot);
        rArmTf.Children.Add(_bodyOff);
        _mdlRightArm!.Transform = rArmTf;

        var lLegTf = new Transform3DGroup();
        lLegTf.Children.Add(new RotateTransform3D(_lLegPitch, lLegPivot));
        lLegTf.Children.Add(bodyRot);
        lLegTf.Children.Add(_bodyOff);
        _mdlLeftLeg!.Transform = lLegTf;

        var rLegTf = new Transform3DGroup();
        rLegTf.Children.Add(new RotateTransform3D(_rLegPitch, rLegPivot));
        rLegTf.Children.Add(bodyRot);
        rLegTf.Children.Add(_bodyOff);
        _mdlRightLeg!.Transform = rLegTf;

        _transformsBuilt = true;
    }

    private void OnAnimTick(object? sender, EventArgs e)
    {
        if (_parts == null) return;
        BuildTransforms();

        double t = (DateTime.Now - _animStart).TotalSeconds;
        var pose = AnimationController.Evaluate(_animState, t);

        // Just update the angle/offset values — no allocations
        _bodyRot.Angle = pose.BodyTilt;
        _bodyOff.OffsetX = pose.BodyOffset.X;
        _bodyOff.OffsetY = pose.BodyOffset.Y;
        _bodyOff.OffsetZ = pose.BodyOffset.Z;
        _headTilt.Angle = pose.HeadTilt;
        _headTurn.Angle = pose.HeadTurn;
        _lArmPitch.Angle = pose.LeftArmPitch;
        _lArmRaise.Angle = -pose.LeftArmRaise;
        _rArmPitch.Angle = pose.RightArmPitch;
        _rArmRaise.Angle = pose.RightArmRaise;
        _lLegPitch.Angle = pose.LeftLegPitch;
        _rLegPitch.Angle = pose.RightLegPitch;
    }

    // ---- Chat ----

    private void UpdateChatSystemPrompt()
    {
        var prompt = ChatService.BuildSystemPrompt(
            _character.Name,
            _character.Class.ToString().ToLowerInvariant(),
            _character.Personality.Creativity,
            _character.Personality.Chattiness,
            _character.Personality.Diligence,
            _character.Directive);
        _chat.SetSystemPrompt(prompt);
    }

    private async void OnChatSend(object sender, RoutedEventArgs e)
    {
        await SendChatMessage();
    }

    private async void OnChatKeyDown(object sender, KeyEventArgs e)
    {
        if (e.Key == Key.Enter)
        {
            e.Handled = true;
            await SendChatMessage();
        }
    }

    private async Task SendChatMessage()
    {
        var text = ChatInput.Text.Trim();
        if (string.IsNullOrEmpty(text)) return;

        ChatInput.Text = "";
        AppendChat($"You: {text}\n", "#88ff88");

        if (!_chat.IsAvailable)
        {
            await CheckChatServerAsync();
            if (!_chat.IsAvailable)
            {
                AppendChat("[Not connected — start llama-server first]\n", "#ff6666");
                return;
            }
        }

        UpdateChatSystemPrompt();

        // Show animated thinking dots
        _isThinking = true;
        _animState = AnimState.Idle;
        _animStart = DateTime.Now;
        var thinkingRun = new System.Windows.Documents.Run($"{_character.Name}: thinking")
        {
            Foreground = new SolidColorBrush(Color.FromRgb(0xff, 0xff, 0x66))
        };
        ChatParagraph.Inlines.Add(thinkingRun);
        ChatHistory.ScrollToEnd();

        // Animate the dots while waiting
        var cts = new CancellationTokenSource();
        _ = AnimateDotsAsync(thinkingRun, _character.Name, cts.Token);

        _chat.SetAgentName(_character.Name);
        var reply = await _chat.SendAsync(text, _character.Name);

        // Stop dots, remove thinking indicator
        cts.Cancel();
        ChatParagraph.Inlines.Remove(thinkingRun);
        _isThinking = false;
        AppendChat($"{_character.Name}: {reply}\n\n", "#cccccc");

        // Play reaction animation, then return to idle after 4 seconds
        _animState = PickAnimForReply(reply);
        _animStart = DateTime.Now;
        _ = ReturnToIdleAfter(4.0);
    }

    private async Task AnimateDotsAsync(System.Windows.Documents.Run run, string name, CancellationToken ct)
    {
        string[] frames = [" .  ", " .. ", " ...", "    "];
        int i = 0;
        while (!ct.IsCancellationRequested)
        {
            try
            {
                await Task.Delay(400, ct);
                Dispatcher.Invoke(() => run.Text = $"{name}: thinking{frames[i++ % frames.Length]}");
            }
            catch { break; }
        }
    }

    private async Task ReturnToIdleAfter(double seconds)
    {
        await Task.Delay(TimeSpan.FromSeconds(seconds));
        if (!_isThinking)
        {
            _animState = AnimState.Idle;
            _animStart = DateTime.Now;
        }
    }

    private static AnimState PickAnimForReply(string reply)
    {
        var lower = reply.ToLowerInvariant();
        if (lower.Contains("hello") || lower.Contains("hey") || lower.Contains("hi ") || lower.Contains("greet"))
            return AnimState.Wave;
        if (lower.Contains("great") || lower.Contains("awesome") || lower.Contains("excellent") ||
            lower.Contains("perfect") || lower.Contains("wonderful") || lower.Contains("!"))
            return AnimState.Celebrate;
        if (lower.Contains("build") || lower.Contains("place") || lower.Contains("construct"))
            return AnimState.Place;
        if (lower.Contains("gather") || lower.Contains("pick") || lower.Contains("collect") || lower.Contains("grab"))
            return AnimState.PickUp;
        if (lower.Contains("go") || lower.Contains("walk") || lower.Contains("head") || lower.Contains("move"))
            return AnimState.Walk;
        return AnimState.Talk;
    }

    private void OnShowAdvanced(object sender, RoutedEventArgs e)
    {
        SyncFromUi();
        var (_, intrinsic, personality, userDir) = ChatService.BuildSystemPromptParts(
            _character.Name,
            _character.Class.ToString().ToLowerInvariant(),
            _character.Personality.Creativity,
            _character.Personality.Chattiness,
            _character.Personality.Diligence,
            _character.Directive);

        var text = $"{intrinsic}\n\n{personality}";
        if (!string.IsNullOrEmpty(userDir))
            text += $"\n\n{userDir}";

        var win = new Window
        {
            Title = "Advanced Directives (Read-Only)",
            Width = 650, Height = 550,
            Background = new SolidColorBrush(Color.FromRgb(0x1a, 0x1a, 0x2e)),
            WindowStartupLocation = WindowStartupLocation.CenterOwner,
            Owner = this,
        };

        var tb = new TextBox
        {
            Text = text,
            IsReadOnly = true,
            AcceptsReturn = true,
            TextWrapping = TextWrapping.Wrap,
            VerticalScrollBarVisibility = ScrollBarVisibility.Auto,
            Background = new SolidColorBrush(Color.FromRgb(0x0a, 0x0a, 0x18)),
            Foreground = new SolidColorBrush(Color.FromRgb(0xcc, 0xcc, 0xcc)),
            BorderThickness = new Thickness(0),
            FontFamily = new FontFamily("Consolas"),
            FontSize = 12,
            Padding = new Thickness(12),
        };
        win.Content = tb;
        win.ShowDialog();
    }

    private void OnChatClear(object sender, RoutedEventArgs e)
    {
        _chat.ClearHistory();
        ChatParagraph.Inlines.Clear();
    }

    private void AppendChat(string text, string colorHex)
    {
        var color = (Color)ColorConverter.ConvertFromString(colorHex);
        var run = new System.Windows.Documents.Run(text)
        {
            Foreground = new SolidColorBrush(color)
        };
        ChatParagraph.Inlines.Add(run);
        ChatHistory.ScrollToEnd();
    }
}
