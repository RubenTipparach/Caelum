using System.Text.Json;
using System.Text.Json.Serialization;

namespace CharEditor.Models;

/// <summary>
/// JSON export format matching what the C game reads from cache/worlds/{id}/agents/
/// </summary>
public class AgentExportData
{
    [JsonPropertyName("id")] public string Id { get; set; } = "";
    [JsonPropertyName("name")] public string Name { get; set; } = "";
    [JsonPropertyName("class")] public string Class { get; set; } = "worker";
    [JsonPropertyName("seed")] public uint Seed { get; set; }
    [JsonPropertyName("appearance")] public AppearanceData Appearance { get; set; } = new();
    [JsonPropertyName("personality")] public PersonalityData Personality { get; set; } = new();
    [JsonPropertyName("attributes")] public AttributeData Attributes { get; set; } = new();
    [JsonPropertyName("directive")] public string Directive { get; set; } = "";

    public class AppearanceData
    {
        [JsonPropertyName("height")] public double Height { get; set; }
        [JsonPropertyName("stockiness")] public double Stockiness { get; set; }
        [JsonPropertyName("head_scale")] public double HeadScale { get; set; }
        [JsonPropertyName("arm_length")] public double ArmLength { get; set; }
        [JsonPropertyName("torso_width")] public double TorsoWidth { get; set; }
        [JsonPropertyName("leg_length")] public double LegLength { get; set; }
        [JsonPropertyName("primary_color")] public double[] PrimaryColor { get; set; } = [0.5, 0.5, 0.6];
        [JsonPropertyName("accent_color")] public double[] AccentColor { get; set; } = [0.8, 0.7, 0.2];
    }

    public class PersonalityData
    {
        [JsonPropertyName("creativity")] public double Creativity { get; set; }
        [JsonPropertyName("chattiness")] public double Chattiness { get; set; }
        [JsonPropertyName("diligence")] public double Diligence { get; set; }
    }

    public class AttributeData
    {
        [JsonPropertyName("strength")] public int Strength { get; set; }
        [JsonPropertyName("perception")] public int Perception { get; set; }
        [JsonPropertyName("endurance")] public int Endurance { get; set; }
    }

    public static AgentExportData FromCharacter(AgentCharacter c)
    {
        return new AgentExportData
        {
            Id = c.Id,
            Name = c.Name,
            Class = c.Class.ToString().ToLowerInvariant(),
            Seed = c.Seed,
            Directive = c.Directive,
            Appearance = new AppearanceData
            {
                Height = Math.Round(c.Appearance.Height, 2),
                Stockiness = Math.Round(c.Appearance.Stockiness, 2),
                HeadScale = Math.Round(c.Appearance.HeadScale, 2),
                ArmLength = Math.Round(c.Appearance.ArmLength, 2),
                TorsoWidth = Math.Round(c.Appearance.TorsoWidth, 2),
                LegLength = Math.Round(c.Appearance.LegLength, 2),
                PrimaryColor = [
                    Math.Round(c.Appearance.PrimaryR, 2),
                    Math.Round(c.Appearance.PrimaryG, 2),
                    Math.Round(c.Appearance.PrimaryB, 2)
                ],
                AccentColor = [
                    Math.Round(c.Appearance.AccentR, 2),
                    Math.Round(c.Appearance.AccentG, 2),
                    Math.Round(c.Appearance.AccentB, 2)
                ]
            },
            Personality = new PersonalityData
            {
                Creativity = Math.Round(c.Personality.Creativity, 2),
                Chattiness = Math.Round(c.Personality.Chattiness, 2),
                Diligence = Math.Round(c.Personality.Diligence, 2)
            },
            Attributes = new AttributeData
            {
                Strength = c.Attributes.Strength,
                Perception = c.Attributes.Perception,
                Endurance = c.Attributes.Endurance
            }
        };
    }

    public static AgentCharacter ToCharacter(AgentExportData d)
    {
        var c = new AgentCharacter
        {
            Id = string.IsNullOrEmpty(d.Id) ? Guid.NewGuid().ToString() : d.Id,
            Name = d.Name,
            Seed = d.Seed,
            Directive = d.Directive
        };
        if (Enum.TryParse<AgentClass>(d.Class, true, out var cls))
            c.Class = cls;

        var a = c.Appearance;
        a.Height = d.Appearance.Height;
        a.Stockiness = d.Appearance.Stockiness;
        a.HeadScale = d.Appearance.HeadScale;
        a.ArmLength = d.Appearance.ArmLength;
        a.TorsoWidth = d.Appearance.TorsoWidth;
        a.LegLength = d.Appearance.LegLength;
        if (d.Appearance.PrimaryColor.Length >= 3)
        {
            a.PrimaryR = d.Appearance.PrimaryColor[0];
            a.PrimaryG = d.Appearance.PrimaryColor[1];
            a.PrimaryB = d.Appearance.PrimaryColor[2];
        }
        if (d.Appearance.AccentColor.Length >= 3)
        {
            a.AccentR = d.Appearance.AccentColor[0];
            a.AccentG = d.Appearance.AccentColor[1];
            a.AccentB = d.Appearance.AccentColor[2];
        }

        c.Personality.Creativity = d.Personality.Creativity;
        c.Personality.Chattiness = d.Personality.Chattiness;
        c.Personality.Diligence = d.Personality.Diligence;
        c.Attributes.Strength = d.Attributes.Strength;
        c.Attributes.Perception = d.Attributes.Perception;
        c.Attributes.Endurance = d.Attributes.Endurance;

        return c;
    }

    public string ToJson()
    {
        return JsonSerializer.Serialize(this, new JsonSerializerOptions
        {
            WriteIndented = true
        });
    }

    public static AgentExportData? FromJson(string json)
    {
        return JsonSerializer.Deserialize<AgentExportData>(json);
    }
}
