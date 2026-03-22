using System.Diagnostics;
using System.IO;
using System.Net.Http;
using System.Text;
using System.Text.Json;

namespace CharEditor.Services;

/// <summary>
/// Communicates with a local llama-server for personality testing.
/// No grammar constraints — free conversation only, no action execution.
/// Can auto-launch llama-server if the binary + model are found.
/// </summary>
public class ChatService
{
    private readonly HttpClient _http;
    private readonly List<ChatMessage> _history = new();
    private string _systemPrompt = "";
    private int _port = 8080;
    private Process? _serverProcess;
    private string? _logPath;

    public record ChatMessage(string Role, string Content);

    public bool IsAvailable { get; private set; }
    public Action<string>? OnStatusUpdate { get; set; }

    public ChatService()
    {
        _http = new HttpClient { Timeout = TimeSpan.FromSeconds(120) };
    }

    public void SetPort(int port)
    {
        _port = port;
    }

    public void SetSystemPrompt(string prompt)
    {
        _systemPrompt = prompt;
    }

    public void ClearHistory()
    {
        _history.Clear();
    }

    /// <summary>
    /// Set the agent name for logging. Creates/appends to cache/agents/{name}/chat_log.jsonl
    /// </summary>
    public void SetAgentName(string name)
    {
        // Find project root
        var dir = AppDomain.CurrentDomain.BaseDirectory;
        for (int i = 0; i < 8; i++)
        {
            var parent = Path.GetDirectoryName(dir);
            if (parent == null) break;
            dir = parent;
            if (File.Exists(Path.Combine(dir, "CMakeLists.txt")))
            {
                var safeName = string.Join("_", name.ToLowerInvariant().Split(Path.GetInvalidFileNameChars()));
                var logDir = Path.Combine(dir, "cache", "agents", safeName);
                Directory.CreateDirectory(logDir);
                _logPath = Path.Combine(logDir, "chat_log.jsonl");
                return;
            }
        }
    }

    private void LogEntry(string type, string agent, string text)
    {
        if (_logPath == null) return;
        try
        {
            var ts = DateTime.UtcNow.ToString("o");
            var escaped = text.Replace("\\", "\\\\").Replace("\"", "\\\"").Replace("\n", "\\n").Replace("\r", "");
            var line = $"{{\"ts\":\"{ts}\",\"type\":\"{type}\",\"agent\":\"{agent}\",\"text\":\"{escaped}\"}}\n";
            File.AppendAllText(_logPath, line);
        }
        catch { /* don't crash chat over logging */ }
    }

    public async Task<bool> CheckHealthAsync()
    {
        try
        {
            var resp = await _http.GetAsync($"http://127.0.0.1:{_port}/health");
            IsAvailable = resp.IsSuccessStatusCode;
            return IsAvailable;
        }
        catch
        {
            IsAvailable = false;
            return false;
        }
    }

    /// <summary>
    /// Try to find and launch llama-server if not already running.
    /// Walks up from exe dir to find the project root (CMakeLists.txt).
    /// </summary>
    public async Task<string> TryAutoLaunchAsync()
    {
        // Already running?
        if (await CheckHealthAsync())
            return "llama-server already running";

        // Find project root
        string? projectRoot = null;
        var dir = AppDomain.CurrentDomain.BaseDirectory;
        for (int i = 0; i < 8; i++)
        {
            var parent = Path.GetDirectoryName(dir);
            if (parent == null) break;
            dir = parent;
            if (File.Exists(Path.Combine(dir, "CMakeLists.txt")))
            {
                projectRoot = dir;
                break;
            }
        }
        if (projectRoot == null)
            return "Could not find project root";

        var serverExe = Path.Combine(projectRoot, "tools", "ai", "llama-server.exe");
        if (!File.Exists(serverExe))
            return $"llama-server not found — run tools\\ai-install.bat";

        // Find model file
        var aiDir = Path.Combine(projectRoot, "tools", "ai");
        var modelFile = Directory.GetFiles(aiDir, "*.gguf").FirstOrDefault();
        if (modelFile == null)
            return "No .gguf model found — run tools\\ai-install.bat";

        // Launch server
        // Launch with visible console window + log stderr to file
        var logFile = Path.Combine(aiDir, "server.log");
        var serverArgs = $"-m \"{modelFile}\" --port {_port} -ngl 20 --ctx-size 16384 --flash-attn on --threads 8 --log-file \"{logFile}\"";
        var psi = new ProcessStartInfo
        {
            FileName = serverExe,
            Arguments = serverArgs,
            WorkingDirectory = aiDir,
            UseShellExecute = true,
        };

        try
        {
            _serverProcess = Process.Start(psi);
            if (_serverProcess == null)
                return "Failed to start llama-server";

            // Wait for server to become healthy (model loading takes time)
            for (int i = 0; i < 120; i++)
            {
                await Task.Delay(1000);

                // Check if process crashed
                if (_serverProcess.HasExited)
                {
                    var exitCode = _serverProcess.ExitCode;
                    // Read last few lines of log
                    var lastLines = "";
                    try { lastLines = File.ReadAllText(logFile); if (lastLines.Length > 500) lastLines = lastLines[^500..]; } catch { }
                    return $"llama-server crashed (exit code {exitCode}). Check {logFile}\n{lastLines}";
                }

                if (await CheckHealthAsync())
                    return $"llama-server started (PID {_serverProcess.Id})";
                OnStatusUpdate?.Invoke($"Loading model... ({i + 1}s)");
            }
            return $"llama-server not responding after 2min. Check {logFile}";
        }
        catch (Exception ex)
        {
            return $"Failed to launch: {ex.Message}";
        }
    }

    public void StopServer()
    {
        if (_serverProcess != null && !_serverProcess.HasExited)
        {
            _serverProcess.Kill();
            _serverProcess.WaitForExit(3000);
            _serverProcess = null;
        }
    }

    public async Task<string> SendAsync(string userMessage, string agentName = "agent")
    {
        _history.Add(new ChatMessage("user", userMessage));
        LogEntry("player_msg", agentName, userMessage);

        // Build messages array
        var messages = new List<object>();
        if (!string.IsNullOrWhiteSpace(_systemPrompt))
            messages.Add(new { role = "system", content = _systemPrompt });
        for (int mi = 0; mi < _history.Count; mi++)
        {
            var msg = _history[mi];
            // Append /no_think to last user message to disable Qwen3 chain-of-thought
            var c = msg.Content;
            if (msg.Role == "user" && mi == _history.Count - 1)
                c += " /no_think";
            messages.Add(new { role = msg.Role, content = c });
        }

        var payload = new
        {
            model = "any",
            messages,
            temperature = 1.05,
            max_tokens = 300,
            stream = false,
            cache_prompt = false,
            seed = -1,
            repeat_penalty = 1.3,
            repeat_last_n = 128
        };

        var json = JsonSerializer.Serialize(payload);
        var content = new StringContent(json, Encoding.UTF8, "application/json");

        try
        {
            var resp = await _http.PostAsync(
                $"http://127.0.0.1:{_port}/v1/chat/completions", content);

            if (!resp.IsSuccessStatusCode)
            {
                var err = await resp.Content.ReadAsStringAsync();
                return $"[Error: HTTP {(int)resp.StatusCode}] {err[..Math.Min(err.Length, 200)]}";
            }

            var body = await resp.Content.ReadAsStringAsync();

            // Parse content from response
            var reply = ExtractContent(body);
            _history.Add(new ChatMessage("assistant", reply));
            LogEntry("agent_msg", agentName, reply);

            // Keep history manageable
            if (_history.Count > 20)
                _history.RemoveRange(0, 2);

            return reply;
        }
        catch (TaskCanceledException)
        {
            return "[Error: Request timed out — is llama-server still loading the model?]";
        }
        catch (HttpRequestException ex)
        {
            return $"[Error: Cannot connect to llama-server on port {_port}. Run ai-install.bat first.]\n{ex.Message}";
        }
    }

    private static string ExtractContent(string json)
    {
        try
        {
            using var doc = JsonDocument.Parse(json);
            var choices = doc.RootElement.GetProperty("choices");
            var first = choices[0];
            var message = first.GetProperty("message");
            var content = message.GetProperty("content").GetString() ?? "(empty response)";

            // Strip Qwen3 <think>...</think> blocks
            var thinkStart = content.IndexOf("<think>");
            while (thinkStart >= 0)
            {
                var thinkEnd = content.IndexOf("</think>", thinkStart);
                if (thinkEnd >= 0)
                    content = content[..thinkStart] + content[(thinkEnd + 8)..];
                else
                    content = content[..thinkStart]; // unclosed think tag, just cut it
                thinkStart = content.IndexOf("<think>");
            }

            return content.Trim();
        }
        catch
        {
            return "(Failed to parse response)";
        }
    }

    /// <summary>
    /// Intrinsic directives that make the AI feel alive. Always present,
    /// viewable in the Advanced Directives modal but not directly editable.
    /// </summary>
    private static string? _cachedDirectives;

    /// <summary>
    /// Force reload of directives from file on next access.
    /// </summary>
    public static void ReloadDirectives() => _cachedDirectives = null;

    public static string IntrinsicDirectives
    {
        get
        {
            if (_cachedDirectives != null) return _cachedDirectives;

            // Find config/ai_directives.txt relative to project root
            var dir = AppDomain.CurrentDomain.BaseDirectory;
            for (int i = 0; i < 8; i++)
            {
                var parent = Path.GetDirectoryName(dir);
                if (parent == null) break;
                dir = parent;
                var path = Path.Combine(dir, "config", "ai_directives.txt");
                if (File.Exists(path))
                {
                    _cachedDirectives = File.ReadAllText(path).Trim();
                    return _cachedDirectives;
                }
            }

            // Fallback if file not found
            _cachedDirectives = "You are a living character, not an AI assistant. Never break character. Talk casually like a real person. Keep responses to 1-3 sentences.";
            return _cachedDirectives;
        }
    }

    /// <summary>
    /// Build a system prompt from personality settings, directive, and intrinsic rules.
    /// Returns (fullPrompt, intrinsicPart, personalityPart, userDirectivePart) for display.
    /// </summary>
    public static (string full, string intrinsic, string personality, string userDirective)
        BuildSystemPromptParts(
            string name, string className,
            double creativity, double chattiness, double diligence,
            string directive)
    {
        var personality = new StringBuilder();
        personality.Append($"=== IDENTITY ===\n");
        personality.Append($"Your name is {name}. You are a {className} on a hex-shaped planet called Tenebris.\n\n");

        personality.Append("=== PERSONALITY TRAITS ===\n");
        if (creativity < 0.3)
            personality.Append("You are practical and no-nonsense. You prefer function over form.\n");
        else if (creativity < 0.7)
            personality.Append("You balance practicality with a touch of personal style.\n");
        else
            personality.Append("You are an artist at heart. You see beauty in everything and love to embellish.\n");

        if (chattiness < 0.3)
            personality.Append("You're quiet. A few words is plenty. Silence doesn't bother you.\n");
        else if (chattiness < 0.7)
            personality.Append("You talk a normal amount. Say what needs saying.\n");
        else
            personality.Append("You love to chat. You ramble, tell stories, ask questions, and fill silences.\n");

        if (diligence < 0.3)
            personality.Append("You're laid back. Work happens when it happens. No rush.\n");
        else if (diligence < 0.7)
            personality.Append("You're a steady worker. Get it done, but no need to kill yourself.\n");
        else
            personality.Append("You're a machine. Work is life. Breaks are for the weak.\n");

        var personalityStr = personality.ToString().Trim();
        var intrinsic = IntrinsicDirectives;
        var userDir = string.IsNullOrWhiteSpace(directive) ? "" : $"=== USER DIRECTIVE ===\n{directive}";

        var full = $"{intrinsic}\n\n{personalityStr}";
        if (!string.IsNullOrEmpty(userDir))
            full += $"\n\n{userDir}";

        return (full, intrinsic, personalityStr, userDir);
    }

    /// <summary>
    /// Convenience: build and return just the full prompt string.
    /// </summary>
    public static string BuildSystemPrompt(
        string name, string className,
        double creativity, double chattiness, double diligence,
        string directive)
    {
        return BuildSystemPromptParts(name, className, creativity, chattiness, diligence, directive).full;
    }
}
