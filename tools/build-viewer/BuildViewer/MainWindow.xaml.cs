using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using System.Windows;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Media3D;
using Microsoft.Win32;

namespace BuildViewer;

public partial class MainWindow : Window
{
    private string _currentFolder = "";
    private Point _lastMouse;
    private bool _rotating;
    private double _camYaw = -45, _camPitch = 30, _camDist = 15;
    private bool _hexMode = true;
    private List<(int q, int r, int layer, string block)> _lastBlocks = new();

    public MainWindow()
    {
        InitializeComponent();

        // Try to auto-open the default orders folder
        var exeDir = AppDomain.CurrentDomain.BaseDirectory;
        for (int i = 0; i < 8; i++)
        {
            var parent = Path.GetDirectoryName(exeDir);
            if (parent == null) break;
            exeDir = parent;
            var ordersDir = Path.Combine(exeDir, "cache", "agents");
            if (Directory.Exists(ordersDir))
            {
                // Find first agent with orders
                foreach (var agentDir in Directory.GetDirectories(ordersDir))
                {
                    var od = Path.Combine(agentDir, "orders");
                    if (Directory.Exists(od))
                    {
                        LoadFolder(od);
                        break;
                    }
                }
                break;
            }
        }

        // Mouse orbit controls — works anywhere on the window
        MouseDown += (s, e) => { _rotating = true; _lastMouse = e.GetPosition(this); CaptureMouse(); };
        MouseUp += (s, e) => { _rotating = false; ReleaseMouseCapture(); };
        MouseMove += (s, e) =>
        {
            if (!_rotating) return;
            var p = e.GetPosition(this);
            _camYaw += (p.X - _lastMouse.X) * 0.5;
            _camPitch += (p.Y - _lastMouse.Y) * 0.5;
            _camPitch = Math.Clamp(_camPitch, -89, 89);
            _lastMouse = p;
            UpdateCamera();
        };
        MouseWheel += (s, e) =>
        {
            _camDist -= e.Delta * 0.01;
            _camDist = Math.Clamp(_camDist, 2, 50);
            UpdateCamera();
        };
    }

    private void UpdateCamera()
    {
        double yr = _camYaw * Math.PI / 180;
        double pr = _camPitch * Math.PI / 180;
        double x = _camDist * Math.Cos(pr) * Math.Sin(yr);
        double y = _camDist * Math.Sin(pr);
        double z = _camDist * Math.Cos(pr) * Math.Cos(yr);
        Camera.Position = new Point3D(x, y, z);
        Camera.LookDirection = new Vector3D(-x, -y, -z);
    }

    private void OnOpenFolder(object sender, RoutedEventArgs e)
    {
        var dlg = new OpenFileDialog
        {
            Title = "Select a build script JSON",
            Filter = "JSON files|*.json",
            InitialDirectory = _currentFolder
        };
        if (dlg.ShowDialog() == true)
        {
            LoadFolder(Path.GetDirectoryName(dlg.FileName)!);
            var name = Path.GetFileNameWithoutExtension(dlg.FileName);
            for (int i = 0; i < ScriptList.Items.Count; i++)
            {
                if ((string)ScriptList.Items[i] == name)
                {
                    ScriptList.SelectedIndex = i;
                    break;
                }
            }
        }
    }

    private void LoadFolder(string folder)
    {
        _currentFolder = folder;
        ScriptList.Items.Clear();
        foreach (var f in Directory.GetFiles(folder, "*.json"))
            ScriptList.Items.Add(Path.GetFileNameWithoutExtension(f));
        Title = $"Build Script Viewer — {folder}";
    }

    private void OnScriptSelected(object sender, System.Windows.Controls.SelectionChangedEventArgs e)
    {
        if (ScriptList.SelectedItem == null) return;
        var name = (string)ScriptList.SelectedItem;
        var path = Path.Combine(_currentFolder, name + ".json");
        if (!File.Exists(path)) return;

        try
        {
            var json = File.ReadAllText(path);
            var doc = JsonDocument.Parse(json);
            var root = doc.RootElement;

            var scriptName = root.TryGetProperty("name", out var n) ? n.GetString() ?? "" : "";
            var desc = root.TryGetProperty("description", out var d) ? d.GetString() ?? "" : "";

            TxtName.Text = scriptName;
            TxtDesc.Text = desc;

            // Parse steps
            var blocks = new List<(int q, int r, int layer, string block)>();
            int moveCount = 0, sayCount = 0, placeCount = 0, breakCount = 0;

            if (root.TryGetProperty("steps", out var steps))
            {
                foreach (var step in steps.EnumerateArray())
                {
                    var action = step.TryGetProperty("action", out var a) ? a.GetString() ?? "" : "";
                    switch (action)
                    {
                        case "place":
                        case "place_rel":
                            int q = step.TryGetProperty("q", out var qv) ? qv.GetInt32() :
                                    step.TryGetProperty("dq", out var dqv) ? dqv.GetInt32() : 0;
                            int r = step.TryGetProperty("r", out var rv) ? rv.GetInt32() :
                                    step.TryGetProperty("dr", out var drv) ? drv.GetInt32() : 0;
                            int l = step.TryGetProperty("layer", out var lv) ? lv.GetInt32() :
                                    step.TryGetProperty("dlayer", out var dlv) ? dlv.GetInt32() : 0;
                            var block = step.TryGetProperty("block", out var bv) ? bv.GetString() ?? "stone" : "stone";
                            blocks.Add((q, r, l, block));
                            placeCount++;
                            break;
                        case "break":
                        case "break_rel":
                            breakCount++;
                            break;
                        case "move_to":
                        case "move_rel":
                            moveCount++;
                            break;
                        case "say":
                            sayCount++;
                            break;
                    }
                }
            }

            TxtStats.Text = $"Steps: {steps.GetArrayLength()}\n" +
                           $"Place: {placeCount}  Break: {breakCount}  Move: {moveCount}  Say: {sayCount}\n" +
                           $"Blocks: {blocks.Count}";

            // Build 3D preview
            BuildPreview(blocks);
        }
        catch (Exception ex)
        {
            TxtName.Text = "Error";
            TxtDesc.Text = ex.Message;
        }
    }

    private static Color BlockColor(string block) => block switch
    {
        "stone" => Color.FromRgb(140, 140, 150),
        "dirt" => Color.FromRgb(139, 90, 43),
        "grass" => Color.FromRgb(76, 153, 0),
        "sand" => Color.FromRgb(210, 190, 130),
        "ice" => Color.FromRgb(160, 210, 240),
        "torch" => Color.FromRgb(255, 180, 50),
        "water" => Color.FromRgb(30, 100, 200),
        _ => Color.FromRgb(180, 180, 180),
    };

    private void BuildPreview(List<(int q, int r, int layer, string block)> blocks)
    {
        _lastBlocks = blocks;

        if (blocks.Count == 0)
        {
            BlocksModel.Content = null;
            return;
        }

        int minQ = blocks.Min(b => b.q), maxQ = blocks.Max(b => b.q);
        int minR = blocks.Min(b => b.r), maxR = blocks.Max(b => b.r);
        int minL = blocks.Min(b => b.layer), maxL = blocks.Max(b => b.layer);
        double cq = (minQ + maxQ) / 2.0;
        double cl = (minL + maxL) / 2.0;
        double cr = (minR + maxR) / 2.0;

        var group = new Model3DGroup();
        double sqrt3 = 1.7320508;

        foreach (var (q, r, layer, block) in blocks)
        {
            var color = BlockColor(block);
            var material = new DiffuseMaterial(new SolidColorBrush(color));

            double hx, hy, hz;
            MeshGeometry3D mesh;

            if (_hexMode)
            {
                // Hex grid: radius=0.9, spacing matches game proportions
                double hexR = 0.9;
                hx = (q - cq) * (1.5 * hexR);
                hy = (layer - cl) * 1.0;
                hz = ((q & 1) != 0 ? (r + 0.5 - cr) : (r - cr)) * (sqrt3 * hexR);
                mesh = MakeHexPrism(hx, hy, hz, hexR, 0.45);
            }
            else
            {
                // Square grid: simple 1:1 spacing
                hx = (q - cq) * 1.0;
                hy = (layer - cl) * 1.0;
                hz = (r - cr) * 1.0;
                mesh = MakeCube(hx, hy, hz, 0.45);
            }

            var model = new GeometryModel3D(mesh, material);
            model.BackMaterial = material;
            group.Children.Add(model);
        }

        BlocksModel.Content = group;

        double size = Math.Max(maxQ - minQ, Math.Max(maxR - minR, maxL - minL));
        _camDist = Math.Max(size * 2.0, 5);
        UpdateCamera();
    }

    private void OnToggleShape(object sender, RoutedEventArgs e)
    {
        _hexMode = !_hexMode;
        BtnToggleShape.Content = _hexMode ? "Shape: Hex" : "Shape: Cube";
        BuildPreview(_lastBlocks);
    }

    private static MeshGeometry3D MakeCube(double cx, double cy, double cz, double half)
    {
        var mesh = new MeshGeometry3D();
        double x0 = cx - half, x1 = cx + half;
        double y0 = cy - half, y1 = cy + half;
        double z0 = cz - half, z1 = cz + half;

        void Face(double ax, double ay, double az, double bx, double by, double bz,
                  double ccx, double ccy, double ccz, double dx, double dy, double dz)
        {
            int i = mesh.Positions.Count;
            mesh.Positions.Add(new Point3D(ax, ay, az));
            mesh.Positions.Add(new Point3D(bx, by, bz));
            mesh.Positions.Add(new Point3D(ccx, ccy, ccz));
            mesh.Positions.Add(new Point3D(dx, dy, dz));
            mesh.TriangleIndices.Add(i); mesh.TriangleIndices.Add(i+1); mesh.TriangleIndices.Add(i+2);
            mesh.TriangleIndices.Add(i); mesh.TriangleIndices.Add(i+2); mesh.TriangleIndices.Add(i+3);
        }

        Face(x0,y0,z0, x1,y0,z0, x1,y1,z0, x0,y1,z0);
        Face(x1,y0,z1, x0,y0,z1, x0,y1,z1, x1,y1,z1);
        Face(x0,y1,z0, x1,y1,z0, x1,y1,z1, x0,y1,z1);
        Face(x0,y0,z1, x1,y0,z1, x1,y0,z0, x0,y0,z0);
        Face(x0,y0,z1, x0,y0,z0, x0,y1,z0, x0,y1,z1);
        Face(x1,y0,z0, x1,y0,z1, x1,y1,z1, x1,y1,z0);

        return mesh;
    }

    // Flat-top hexagonal prism matching the game's hex voxels
    private static MeshGeometry3D MakeHexPrism(double cx, double cy, double cz, double radius, double halfH)
    {
        var mesh = new MeshGeometry3D();

        // 6 vertices for flat-top hexagon
        var hex = new Point[6];
        for (int i = 0; i < 6; i++)
        {
            double angle = Math.PI / 3.0 * i; // flat-top: starts at 0
            hex[i] = new Point(cx + radius * Math.Cos(angle), cz + radius * Math.Sin(angle));
        }

        double yBot = cy - halfH;
        double yTop = cy + halfH;

        // Top face (fan from center)
        int baseIdx = mesh.Positions.Count;
        mesh.Positions.Add(new Point3D(cx, yTop, cz)); // center
        for (int i = 0; i < 6; i++)
            mesh.Positions.Add(new Point3D(hex[i].X, yTop, hex[i].Y));
        for (int i = 0; i < 6; i++)
        {
            mesh.TriangleIndices.Add(baseIdx);
            mesh.TriangleIndices.Add(baseIdx + 1 + i);
            mesh.TriangleIndices.Add(baseIdx + 1 + (i + 1) % 6);
        }

        // Bottom face (fan, reversed winding)
        baseIdx = mesh.Positions.Count;
        mesh.Positions.Add(new Point3D(cx, yBot, cz)); // center
        for (int i = 0; i < 6; i++)
            mesh.Positions.Add(new Point3D(hex[i].X, yBot, hex[i].Y));
        for (int i = 0; i < 6; i++)
        {
            mesh.TriangleIndices.Add(baseIdx);
            mesh.TriangleIndices.Add(baseIdx + 1 + (i + 1) % 6);
            mesh.TriangleIndices.Add(baseIdx + 1 + i);
        }

        // 6 side faces (2 triangles each)
        for (int i = 0; i < 6; i++)
        {
            int j = (i + 1) % 6;
            baseIdx = mesh.Positions.Count;
            mesh.Positions.Add(new Point3D(hex[i].X, yBot, hex[i].Y));
            mesh.Positions.Add(new Point3D(hex[j].X, yBot, hex[j].Y));
            mesh.Positions.Add(new Point3D(hex[j].X, yTop, hex[j].Y));
            mesh.Positions.Add(new Point3D(hex[i].X, yTop, hex[i].Y));
            mesh.TriangleIndices.Add(baseIdx);
            mesh.TriangleIndices.Add(baseIdx + 1);
            mesh.TriangleIndices.Add(baseIdx + 2);
            mesh.TriangleIndices.Add(baseIdx);
            mesh.TriangleIndices.Add(baseIdx + 2);
            mesh.TriangleIndices.Add(baseIdx + 3);
        }

        return mesh;
    }
}
