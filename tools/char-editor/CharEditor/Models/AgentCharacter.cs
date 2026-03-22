using System.ComponentModel;
using System.Text.Json.Serialization;

namespace CharEditor.Models;

public class AgentAppearance : INotifyPropertyChanged
{
    private double _height = 1.0;
    private double _stockiness = 0.5;
    private double _headScale = 1.0;
    private double _armLength = 1.0;
    private double _torsoWidth = 1.0;
    private double _legLength = 1.0;
    private double _primaryR = 0.5, _primaryG = 0.5, _primaryB = 0.6;
    private double _accentR = 0.8, _accentG = 0.7, _accentB = 0.2;

    public double Height { get => _height; set { _height = value; OnChanged(); } }
    public double Stockiness { get => _stockiness; set { _stockiness = value; OnChanged(); } }
    public double HeadScale { get => _headScale; set { _headScale = value; OnChanged(); } }
    public double ArmLength { get => _armLength; set { _armLength = value; OnChanged(); } }
    public double TorsoWidth { get => _torsoWidth; set { _torsoWidth = value; OnChanged(); } }
    public double LegLength { get => _legLength; set { _legLength = value; OnChanged(); } }

    // Primary color (body)
    public double PrimaryR { get => _primaryR; set { _primaryR = value; OnChanged(); } }
    public double PrimaryG { get => _primaryG; set { _primaryG = value; OnChanged(); } }
    public double PrimaryB { get => _primaryB; set { _primaryB = value; OnChanged(); } }

    // Accent color (trim/detail)
    public double AccentR { get => _accentR; set { _accentR = value; OnChanged(); } }
    public double AccentG { get => _accentG; set { _accentG = value; OnChanged(); } }
    public double AccentB { get => _accentB; set { _accentB = value; OnChanged(); } }

    public event PropertyChangedEventHandler? PropertyChanged;
    private void OnChanged([System.Runtime.CompilerServices.CallerMemberName] string? name = null)
        => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
}

public class AgentPersonality : INotifyPropertyChanged
{
    private double _creativity = 0.5;
    private double _chattiness = 0.5;
    private double _diligence = 0.7;

    public double Creativity { get => _creativity; set { _creativity = value; OnChanged(); } }
    public double Chattiness { get => _chattiness; set { _chattiness = value; OnChanged(); } }
    public double Diligence { get => _diligence; set { _diligence = value; OnChanged(); } }

    public event PropertyChangedEventHandler? PropertyChanged;
    private void OnChanged([System.Runtime.CompilerServices.CallerMemberName] string? name = null)
        => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
}

public class AgentAttributes : INotifyPropertyChanged
{
    private int _strength = 5;
    private int _perception = 5;
    private int _endurance = 5;

    public int Strength { get => _strength; set { _strength = value; OnChanged(); } }
    public int Perception { get => _perception; set { _perception = value; OnChanged(); } }
    public int Endurance { get => _endurance; set { _endurance = value; OnChanged(); } }

    public event PropertyChangedEventHandler? PropertyChanged;
    private void OnChanged([System.Runtime.CompilerServices.CallerMemberName] string? name = null)
        => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
}

public enum AgentClass
{
    Worker,
    Wanderer,
    Villager
}

public class AgentCharacter : INotifyPropertyChanged
{
    private string _id = Guid.NewGuid().ToString();
    private string _name = "Worker";
    private AgentClass _class = AgentClass.Worker;
    private uint _seed;
    private string _directive = "";

    public string Id { get => _id; set { _id = value; OnChanged(); } }
    public string Name { get => _name; set { _name = value; OnChanged(); } }
    public AgentClass Class { get => _class; set { _class = value; OnChanged(); } }
    public uint Seed { get => _seed; set { _seed = value; OnChanged(); } }
    public string Directive { get => _directive; set { _directive = value; OnChanged(); } }

    public AgentAppearance Appearance { get; set; } = new();
    public AgentPersonality Personality { get; set; } = new();
    public AgentAttributes Attributes { get; set; } = new();

    public event PropertyChangedEventHandler? PropertyChanged;
    private void OnChanged([System.Runtime.CompilerServices.CallerMemberName] string? name = null)
        => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));

    public void Randomize(Random rng)
    {
        Seed = (uint)rng.Next();
        Appearance.Height = 0.8 + rng.NextDouble() * 0.4;
        Appearance.Stockiness = rng.NextDouble();
        Appearance.HeadScale = 0.8 + rng.NextDouble() * 0.4;
        Appearance.ArmLength = 0.8 + rng.NextDouble() * 0.4;
        Appearance.TorsoWidth = 0.8 + rng.NextDouble() * 0.4;
        Appearance.LegLength = 0.8 + rng.NextDouble() * 0.4;
        Appearance.PrimaryR = rng.NextDouble();
        Appearance.PrimaryG = rng.NextDouble();
        Appearance.PrimaryB = rng.NextDouble();
        Appearance.AccentR = rng.NextDouble();
        Appearance.AccentG = rng.NextDouble();
        Appearance.AccentB = rng.NextDouble();
        Personality.Creativity = rng.NextDouble();
        Personality.Chattiness = rng.NextDouble();
        Personality.Diligence = rng.NextDouble();
        Attributes.Strength = rng.Next(1, 11);
        Attributes.Perception = rng.Next(1, 11);
        Attributes.Endurance = rng.Next(1, 11);
    }
}
