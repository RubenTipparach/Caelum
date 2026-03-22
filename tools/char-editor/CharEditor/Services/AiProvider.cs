using System.IO;
using System.Net.Http;
using System.Text;
using System.Text.Json;

namespace CharEditor.Services;

public enum AiProviderType
{
    Local,
    Claude
}

/// <summary>
/// Reads config/ai_config.yaml to determine provider settings.
/// </summary>
public static class AiConfig
{
    public static AiProviderType Provider { get; private set; } = AiProviderType.Local;
    public static string ClaudeApiKey { get; private set; } = "";
    public static string ClaudeModel { get; private set; } = "claude-sonnet-4-20250514";
    public static int ClaudeMaxTokens { get; private set; } = 300;
    public static double ClaudeTemperature { get; private set; } = 1.0;

    public static void Load()
    {
        var dir = AppDomain.CurrentDomain.BaseDirectory;
        for (int i = 0; i < 8; i++)
        {
            var parent = Path.GetDirectoryName(dir);
            if (parent == null) break;
            dir = parent;
            var path = Path.Combine(dir, "config", "ai_config.yaml");
            if (File.Exists(path))
            {
                ParseYaml(File.ReadAllText(path));
                return;
            }
        }
    }

    private static void ParseYaml(string yaml)
    {
        foreach (var rawLine in yaml.Split('\n'))
        {
            var line = rawLine.Trim();
            if (line.StartsWith('#') || !line.Contains(':')) continue;

            var parts = line.Split(':', 2);
            var key = parts[0].Trim();
            var val = parts[1].Trim().Trim('"');

            switch (key)
            {
                case "provider":
                    Provider = val.ToLowerInvariant() == "claude"
                        ? AiProviderType.Claude : AiProviderType.Local;
                    break;
                case "api_key": ClaudeApiKey = val; break;
                case "model" when Provider == AiProviderType.Claude: ClaudeModel = val; break;
                case "max_tokens" when Provider == AiProviderType.Claude:
                    if (int.TryParse(val, out var mt)) ClaudeMaxTokens = mt;
                    break;
                case "temperature" when Provider == AiProviderType.Claude:
                    if (double.TryParse(val, out var ct)) ClaudeTemperature = ct;
                    break;
            }
        }
    }
}

/// <summary>
/// Claude API provider — drop-in replacement for local llama-server.
/// </summary>
public class ClaudeProvider
{
    private readonly HttpClient _http;
    private readonly List<ChatService.ChatMessage> _history = new();
    private string _systemPrompt = "";

    public ClaudeProvider()
    {
        _http = new HttpClient
        {
            Timeout = TimeSpan.FromSeconds(60),
            BaseAddress = new Uri("https://api.anthropic.com")
        };
    }

    public void SetSystemPrompt(string prompt) => _systemPrompt = prompt;

    public void ClearHistory() => _history.Clear();

    public async Task<string> SendAsync(string userMessage)
    {
        if (string.IsNullOrEmpty(AiConfig.ClaudeApiKey))
            return "[Error: Claude API key not set in config/ai_config.yaml]";

        _history.Add(new ChatService.ChatMessage("user", userMessage));

        var messages = new List<object>();
        foreach (var msg in _history)
            messages.Add(new { role = msg.Role, content = msg.Content });

        var payload = new
        {
            model = AiConfig.ClaudeModel,
            max_tokens = AiConfig.ClaudeMaxTokens,
            temperature = AiConfig.ClaudeTemperature,
            system = _systemPrompt,
            messages
        };

        var json = JsonSerializer.Serialize(payload);
        var request = new HttpRequestMessage(HttpMethod.Post, "/v1/messages")
        {
            Content = new StringContent(json, Encoding.UTF8, "application/json")
        };
        request.Headers.Add("x-api-key", AiConfig.ClaudeApiKey);
        request.Headers.Add("anthropic-version", "2023-06-01");

        try
        {
            var resp = await _http.SendAsync(request);
            var body = await resp.Content.ReadAsStringAsync();

            if (!resp.IsSuccessStatusCode)
                return $"[Error: HTTP {(int)resp.StatusCode}] {body[..Math.Min(body.Length, 200)]}";

            var reply = ExtractClaudeContent(body);
            _history.Add(new ChatService.ChatMessage("assistant", reply));

            if (_history.Count > 20)
                _history.RemoveRange(0, 2);

            return reply;
        }
        catch (TaskCanceledException)
        {
            return "[Error: Claude API request timed out]";
        }
        catch (HttpRequestException ex)
        {
            return $"[Error: Cannot reach Claude API] {ex.Message}";
        }
    }

    private static string ExtractClaudeContent(string json)
    {
        try
        {
            using var doc = JsonDocument.Parse(json);
            var content = doc.RootElement.GetProperty("content");
            var first = content[0];
            return first.GetProperty("text").GetString() ?? "(empty response)";
        }
        catch
        {
            return "(Failed to parse Claude response)";
        }
    }
}
