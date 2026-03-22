# AI NPC Architecture — Local Real-Time LLM Integration

## Vision

A Minecraft-style NPC bot the player walks up to and has full engineering conversations with. It remembers past conversations, plans structures, and executes block place/break actions autonomously. Desktop-only for now (local llama.cpp); web/mobile AI is disabled (future cloud proxy approach TBD).

---

## Inference Engine Comparison

| Engine | Language | C API | Desktop | Browser | Mobile | License | Verdict |
|--------|----------|-------|---------|---------|--------|---------|---------|
| **llama.cpp** | C/C++ | **Native** | Best | Slow WASM | Partial | MIT | **Primary choice** |
| **WebLLM (MLC)** | JS/WebGPU | No (JS) | N/A | **Best** | N/A | Apache 2.0 | **Web choice** |
| **Ollama** | Go wrapper | HTTP only | Good | No | No | MIT | Dev/prototyping only |
| **whisper.cpp** | C/C++ | **Native** | Best | WASM ok | Ok | MIT | **Voice input** |
| ONNX Runtime | C/C++ | Yes | Ok | WebGPU | Yes | MIT | Backup web path |
| ExecuTorch | C++ | C++ | N/A | N/A | **Best** | BSD | Native mobile only |
| MLX | Python/Swift | No | Mac only | No | No | Apache 2.0 | Skip (no C API) |
| TensorRT-LLM | Python | No | NVIDIA only | No | No | Apache 2.0 | Skip (overkill) |
| PowerInfer | C++ | Fork | Large models | No | No | Apache 2.0 | Skip (for 13B+) |

**Winner: llama.cpp for desktop, WebLLM for browser, whisper.cpp for voice.**

---

## Small Model Comparison (Q4 Quantized)

| Model | Params | Q4 Size | RAM | Tool Use Quality | tok/s (GPU) | License | Best For |
|-------|--------|---------|-----|-----------------|-------------|---------|----------|
| **Qwen3-4B** | 4B | 2.8GB | ~4.5GB | Excellent | 80-150 | **Apache 2.0** | Overall best |
| **Ministral-3B** | 3B | 2.0GB | ~3.5GB | **Best-in-class** | 80-150 | Apache 2.0 | Reliable tool calls |
| **Phi-4-mini** | 3.8B | 2.6GB | ~4GB | Great | 60-120 | MIT | Best reasoning |
| Gemma 3 4B QAT | 4B | 2.6GB | ~4GB | Good | 60-120 | Gemma License | Code/JSON output |
| SmolLM3-3B | 3B | 2.0GB | ~3.5GB | Good | 80-150 | Apache 2.0 | Conversation quality |
| Llama 3.2 3B | 3B | 2.0GB | ~3.5GB | 67% BFCL | 80-150 | Llama License | Safe default |
| Llama 3.2 1B | 1B | 0.7GB | ~1.5GB | Weak | 100-200+ | Llama License | Ultra-light fallback |
| SmolLM2-1.7B | 1.7B | 1.1GB | ~2GB | Limited | 100-200 | Apache 2.0 | Web/mobile chat |
| Gemma 3 1B | 1B | 0.5GB | ~1GB | Basic | 100-200+ | Gemma License | Minimum viable web |
| Qwen3-0.6B | 0.6B | 0.4GB | ~1GB | Minimal | 150-300 | Apache 2.0 | Absolute minimum |

---

## Recommended Model Per Tier

| Tier | Model | Size | Why |
|------|-------|------|-----|
| **Desktop (full)** | Qwen3-4B or Ministral-3B | 2-2.8GB | Best reasoning + tool use, Apache 2.0 |
| **Desktop (lite)** | Llama 3.2 3B | 2GB | Good fallback for weaker GPUs |
| **Web browser** | SmolLM2-1.7B or Gemma 3 1B | 0.5-1.1GB | Fits browser memory budget |
| **Mobile** | Llama 3.2 1B | 0.7GB | Meta-optimized for mobile silicon |

---

## Top Model Deep Dive

Head-to-head comparison of the top contenders for the NPC bot role:

| | **Qwen3-4B** | **Ministral-3B** | **Phi-4-mini (3.8B)** | **DeepSeek-R1-Distill (1.5B)** | **SmolLM3-3B** |
|---|---|---|---|---|---|
| **Params** | 4B | 3B | 3.8B | 1.5B | 3B |
| **Q4 Size** | 2.8GB | 2.0GB | 2.6GB | ~1GB | 2.0GB |
| **RAM** | ~4.5GB | ~3.5GB | ~4GB | ~2GB | ~3.5GB |
| **License** | Apache 2.0 | Apache 2.0 | MIT | MIT | Apache 2.0 |
| **Tool/Function Calling** | Excellent | **Best-in-class** | Great | Weak | Good |
| **Reasoning** | **Best** (rivals 72B) | Good | Very strong | Strong (distilled from R1) | Good |
| **Conversation Quality** | Great | Good | Good (can feel dry) | Poor (feels "thinky") | **Best** |
| **Dual-mode** | Yes (`/think` + `/no_think`) | No | No | Always chain-of-thought | Yes |
| **llama.cpp support** | Full, official GGUF | Full | Verify (was pending) | Full | Full |
| **GBNF grammar compat** | Yes | Yes | Yes | Yes | Yes |

### What Each Model Excels At

**Qwen3-4B** — The all-rounder. Dual-mode is killer for our NPC use case: `/no_think` for fast casual chat ("Hey, how's it going?"), `/think` for complex planning ("Build me a 3-story tower with a spiral staircase"). Apache 2.0 means zero licensing headaches. Alibaba actively maintains it.

**Ministral-3B** — Purpose-built for tool calling. If the main concern is "will the NPC reliably output valid JSON actions without hallucinating garbage," this is the safest bet. Smaller than Qwen3-4B too (2GB vs 2.8GB). Slightly weaker at freeform conversation.

**Phi-4-mini** — Microsoft's math/reasoning beast. Scores like a 7-9B model on reasoning benchmarks. Great if the NPC needs to do complex spatial calculations. Downside: llama.cpp support was being finalized — needs verification. Conversation can feel a bit robotic.

**DeepSeek-R1-Distill 1.5B** — Interesting niche. Tiny (1GB) and punches above its weight on multi-step reasoning because it was distilled from DeepSeek's massive R1 reasoning model. But it always "thinks out loud" (chain-of-thought), making it slow and awkward for casual chat. Best used as a **background planning model**, not the conversation face. Could power a dual-model architecture where one model plans and another talks.

**SmolLM3-3B** — Best pure conversationalist. If the NPC needs to feel natural and human-like in dialogue, this wins. Tool calling is decent but not as reliable as Ministral or Qwen. HuggingFace-backed, Apache 2.0.

### Model Recommendation

**Primary: Qwen3-4B** — the dual-mode (`/think` vs `/no_think`) maps perfectly to the NPC's two jobs:
- Player says "Hey what's up?" -> `/no_think` -> instant casual response
- Player says "Build a defensive wall around this base" -> `/think` -> chain-of-thought planning -> structured action sequence

**Fallback for weaker hardware: Ministral-3B** — 800MB smaller, best tool-call reliability, good enough conversation.

**Optional dual-model architecture:** Use DeepSeek-R1-Distill 1.5B as a background "planning brain" that generates build plans, while Qwen3-4B or SmolLM3-3B handles the player-facing dialogue. Adds complexity but separates concerns cleanly.

---

## Agent Behavior & Physical Presence

### Character Classes

| Class | Role | Autonomy | First Priority |
|-------|------|----------|----------------|
| **Worker** | Personal engineering assistant | Low — obeys player commands | Phase 1 (implement first) |
| **Wanderer** | Roams the world, trades with NPCs | Medium — follows trade routes, interacts independently | Phase 3+ |
| **Villager** | Builds their own settlements | High — self-directed building, forms communities | Phase 3+ |

Workers are the first class to implement. They generally obey commands but have their own personality and can make suggestions.

### Physical Rules

Agents obey the same physical laws as the player:
- **No teleportation.** Agents must walk/pathfind to destinations using A* on the hex grid.
- **Sleep when player is far away.** If the player moves beyond a configurable range (e.g. 200m), the worker enters sleep state. Actions pause, state is saved. This is a hard design constraint until a background simulation system is built.
- **Collision.** Agents use the same swept AABB collision as the player against the voxel grid.
- **Gravity.** Agents fall, orient to the local surface normal, and respect SOI transitions (moon/planet gravity).

### Procedural Character Model

Each agent has a procedurally generated character mesh, cached on first spawn:
- **Body:** Simple humanoid built from hex prism primitives (torso, limbs, head) — same voxel aesthetic as the world
- **Variation:** Height, proportions, color palette seeded from agent ID
- **Personality visualization:** Subtle visual cues tied to personality (e.g. stocky + muted colors = pragmatic worker, lanky + bright = creative)
- **Rendering:** Uses existing `LodVertex` format (pos, normal, color) — no new pipeline needed, same as moon/planet meshes
- **Animation:** Simple procedural — bobbing while walking, arm swing, head turn toward player when in conversation range

API sketch for character generation:
```c
typedef struct {
    uint32_t seed;           // deterministic generation from seed
    float    height;         // 0.8 - 1.2 (relative scale)
    float    stockiness;     // 0.0 (lean) - 1.0 (broad)
    HMM_Vec3 primary_color;  // body color
    HMM_Vec3 accent_color;   // trim/detail color
    AgentClass class_type;   // AGENT_WORKER, AGENT_WANDERER, AGENT_VILLAGER
} AgentAppearance;

typedef struct {
    float creativity;        // 0.0 (strict to spec) - 1.0 (embellishes freely)
    float chattiness;        // how often they speak unprompted
    float diligence;         // work speed vs break frequency
    char  name[32];          // generated or player-assigned
} AgentPersonality;

// Generate and cache the character mesh
sg_buffer agent_mesh_generate(const AgentAppearance* appearance);

// Optional: external C# tool for visual preview/tweaking
// (separate app, exports AgentAppearance as JSON for import into game)
```

Future: a small C# WPF/WinForms app for visually tweaking `AgentAppearance` + `AgentPersonality` with a 3D preview, then exporting as JSON that the game imports. Not needed for Phase 1.

### Creativity Slider

Workers have a `creativity` float (0.0 - 1.0) that affects the LLM system prompt:
- **0.0 (Pragmatic):** "Build exactly to specification. Use minimal materials. No decoration."
- **0.5 (Balanced):** "Follow the general plan but add minor structural improvements where sensible."
- **1.0 (Creative):** "Interpret the request freely. Add decorative elements, vary materials, embellish the design."

This maps directly to the LLM temperature and system prompt wording. At 0.0 creativity, the model runs with low temperature and strict instructions. At 1.0, higher temperature and an open-ended system prompt.

### Agent Vision — World Perception

Agents perceive the world through **structured geometry data**, not rendered images. Image-based vision would require multimodal models (much larger, slower at small sizes) and waste GPU on rendering frames. Text-based spatial data works perfectly with text-only LLMs.

#### Distance-Based Detail Tiers

| Range | Data Fed to Agent | Source |
|-------|-------------------|--------|
| **0-10m** | Full cell detail: hex positions, heights, block types, neighbors, player-placed vs natural | `hex_terrain` chunk data |
| **10-50m** | Summarized regions: "hill to the north (5 cells, avg height 8), flat plains east, water body south" | Aggregated from chunks |
| **50-100m** | Terrain profile: "elevated ridge NW rising to ~15m, valley SE dropping to ~2m" | Simplified chunk stats |
| **100m+** | Coarse elevation data from LOD mesh vertices | Already computed by LOD system |

The 0-10m tier is a serialized hex grid snapshot:
```json
{
  "nearby_cells": [
    { "q": 0, "r": 1, "height": 5, "type": "stone", "placed_by": "player" },
    { "q": 1, "r": 0, "height": 3, "type": "grass", "placed_by": "natural" }
  ],
  "summary_10_50m": "Hill N (avg h=8), flat E, water S",
  "summary_50_100m": "Ridge NW, valley SE",
  "agent_pos": { "q": 0, "r": 0, "height": 4 }
}
```

This is cheap to generate (iterate nearby chunks, already loaded), fits easily in the LLM context window, and gives the agent enough spatial awareness to plan builds and navigate.

---

## Tool Use & Structured Output

### Action Vocabulary

The NPC has a small, well-defined set of actions:

```
place_block(x, y, z, type)
break_block(x, y, z)
say("message")
move_to(x, y, z)
look_at(x, y, z)
```

### GBNF Grammar-Constrained Decoding

llama.cpp supports GBNF grammars that force the model to only output tokens matching valid actions. Invalid JSON is literally impossible at the token-sampling level. Even a 3B model achieves near-100% structured output reliability with this approach.

XGrammar (integrated into llama.cpp 2025) provides up to 100x speedup over naive grammar constraint methods, making this negligible overhead.

### Function Calling Reliability (BFCL Rankings)

1. **Ministral-3B** — purpose-built for native function calling
2. **Qwen3-4B** — Hermes-style tool calling, dual-mode reasoning
3. **Phi-4-mini** — native tool tokens (`<|tool|>` / `</|tool|>`)
4. **Gemma 3 4B** — decent via system prompt

For a small, well-defined tool vocabulary like ours, grammar constraints + any 3B+ model = near-perfect reliability.

---

## System Architecture

### Desktop

```
+-----------------------------------------------------+
|                   hex-planets.exe                    |
|                                                      |
|  +----------+  +-----------+  +------------------+  |
|  | Game Loop |  | AI Client |  | Action Executor  |  |
|  | (render,  |--| (WinHTTP  |--| (place/break/    |  |
|  |  input)   |  |  to local |  |  move/say queue) |  |
|  +----------+  |  server)  |  +------------------+  |
|                 +-----+-----+                        |
+-------------------+---+------------------------------+
                    |
                    | HTTP localhost:8080
                    |
+-------------------+----------------------------------+
|              llama-server (sidecar)                   |
|  +------------------------------------------------+  |
|  |  Qwen3-4B Q4_K_M  |  GBNF grammar enforced    |  |
|  |  /v1/chat/completions (OpenAI-compatible)      |  |
|  +------------------------------------------------+  |
+------------------------------------------------------+

+------------------------------------------------------+
|              whisper-server (optional)                |
|  whisper-tiny/base | microphone -> text              |
+------------------------------------------------------+
```

### Desktop Flow

1. `ai-install.bat` downloads llama.cpp release + GGUF model (~3GB one-time download)
2. Game launches `llama-server` as a child process on startup
3. Player approaches NPC -> chat UI opens
4. Player text (or whisper.cpp voice) -> sent to llama-server via WinHTTP
5. Response parsed as JSON action + dialogue
6. Action executor queues block operations
7. NPC speaks dialogue in chat bubble

### Web/Mobile

AI is **disabled** for web and mobile builds. The NPC entity exists but uses scripted dialogue trees and fixed build recipes instead of LLM inference. A cloud proxy approach (Vercel endpoint -> hosted LLM API) may be explored later as an opt-in feature.

---

## Memory System

### Short-Term Memory
- Rolling conversation in the LLM context window (last N turns)
- Managed via llama.cpp KV cache reuse across turns

### Long-Term Memory
- SQLite database (already used for world persistence)
- Store conversation summaries, NPC state, building project progress
- On new conversation: load relevant past context from SQLite into system prompt

### Action Log
- Track what the NPC has built so it can reference and continue projects
- Stored per-world in the existing save system (`cache/worlds/{id}/`)

---

## NPC Response Format

Each response from the LLM is a JSON object constrained by GBNF grammar:

```json
{
  "dialogue": "I'll build a watchtower here. Starting with the foundation.",
  "actions": [
    { "type": "place_block", "x": 10, "y": 0, "z": 5, "block": "stone" },
    { "type": "place_block", "x": 11, "y": 0, "z": 5, "block": "stone" },
    { "type": "place_block", "x": 10, "y": 1, "z": 5, "block": "stone" },
    { "type": "say", "text": "Foundation laid. Want me to continue upward?" }
  ],
  "memory_note": "Started watchtower project at (10,0,5). Foundation is 2x1 stone."
}
```

---

## New Files

```
src/ai_npc.h/c          -- AI client (WinHTTP calls to llama-server)
src/ai_actions.h/c       -- Action parser + executor (JSON -> game commands)
src/ai_memory.h/c        -- SQLite-backed conversation memory
src/ai_agent.h/c         -- Agent entity (appearance, personality, state machine, physics)
src/ai_vision.h/c        -- World perception (geometry serialization by distance tier)
src/ai_pathfind.h/c      -- A* pathfinding on hex grid
tools/ai-install.bat     -- Downloads llama-server + model
tools/ai-grammar.gbnf    -- GBNF grammar for valid NPC actions
```

---

## Implementation Outline

### Phase 1: Wire Up AI Model

Get the LLM running and talking to the game.

- Write `ai-install.bat` to download llama-server + GGUF model (Qwen3-4B Q4_K_M, ~3GB)
- Implement `ai_npc.c` — WinHTTP client to localhost llama-server (`/v1/chat/completions`)
- Define GBNF grammar for action vocabulary (place/break/say/move/look)
- Launch llama-server as child process on game startup, kill on shutdown
- Basic in-game chat UI: text input field, chat bubble output above NPC
- Test with Ollama during development for fast prompt iteration
- Verify structured JSON output works reliably with grammar constraints

**Deliverable:** Player can type a message, get a response from the local LLM, see it in-game.

### Phase 2: Character Creator (C# App)

Separate desktop tool for designing and configuring AI characters.

- C# WPF/WinForms app with 3D preview (could use MonoGame, Helix Toolkit, or simple OpenGL control)
- Procedural character mesh preview — hex prism humanoid built from voxel-style primitives
- **Appearance controls:**
  - Height, stockiness, proportions sliders
  - Primary/accent color pickers
  - Body part scaling (head size, arm length, torso width)
  - Randomize button (seed-based, deterministic)
- **Personality/Directive controls:**
  - Name field (editable)
  - Character class dropdown (Worker, Wanderer, Villager)
  - Creativity slider (0.0 pragmatic - 1.0 creative)
  - Chattiness slider (how often they speak unprompted)
  - Diligence slider (work speed vs break frequency)
  - Custom directive text box — free-form instructions injected into the LLM system prompt
    (e.g. "You are a stonemason who prefers granite. You speak in short sentences.")
  - Template directives — presets the user can start from and customize:
    - "Obedient Builder" — follows orders exactly, minimal chatter
    - "Creative Architect" — interprets freely, suggests improvements
    - "Chatty Companion" — talks a lot, asks questions, slower to work
    - "Silent Workhorse" — rarely speaks, high diligence, low creativity
- **Attribute system:**
  - Strength (carry capacity, build speed)
  - Perception (vision range, detail level)
  - Endurance (how long before needing rest)
  - Custom attributes (extensible, key-value pairs for future use)
- **Export:** saves character as JSON file
  ```json
  {
    "name": "Grok",
    "class": "worker",
    "seed": 48291,
    "appearance": {
      "height": 1.05,
      "stockiness": 0.7,
      "primary_color": [0.6, 0.4, 0.3],
      "accent_color": [0.8, 0.7, 0.2]
    },
    "personality": {
      "creativity": 0.3,
      "chattiness": 0.2,
      "diligence": 0.9
    },
    "attributes": {
      "strength": 8,
      "perception": 5,
      "endurance": 7
    },
    "directive": "You are a stonemason. Build with stone whenever possible. Keep responses brief."
  }
  ```
- **Import in game:** C game reads JSON from `cache/worlds/{id}/agents/`, generates mesh from appearance data, injects directive into LLM system prompt

**Deliverable:** Standalone app where you design a character, tweak its personality/directives, preview it in 3D, export JSON that the game imports.

### Phase 3: Game Integration + Interaction Logging

Bring the character into the game world and log everything for analysis.

#### 3a: Physical Agent in Game
- Implement `ai_agent.c` — agent entity with state machine (idle, walking, working, sleeping, conversing)
- Load character JSON on world load, generate procedural mesh from appearance data (cached)
- Agent renders using existing `LodVertex` pipeline (pos, normal, color)
- A* pathfinding on hex grid (`ai_pathfind.c`) — no teleportation
- Collision using existing swept AABB system
- Gravity + SOI awareness (walk on moons, orient to surface normal)
- Sleep when player is beyond 200m (hard constraint until background sim is built)
- Basic procedural animation (walk bob, arm swing, head turn toward player in conversation range)
- Creativity slider maps to LLM temperature + system prompt wording

#### 3b: World Perception (Vision)
- Implement `ai_vision.c` — serialize nearby geometry as structured text
- Distance-based detail tiers fed into LLM context per request:
  - 0-10m: full cell data (positions, heights, block types)
  - 10-50m: summarized regions
  - 50-100m: terrain profile
  - 100m+: LOD mesh vertex data
- Agent references terrain in conversation ("I see a ridge to the north")

#### 3c: Interaction Logging
All agent activity is logged for analysis. Two log types:

**Player-Agent Interaction Log** (`cache/worlds/{id}/ai/chat_log.jsonl`):
```jsonl
{"ts":"2026-03-21T14:30:00Z","type":"player_msg","agent":"Grok","text":"Build a wall here"}
{"ts":"2026-03-21T14:30:01Z","type":"agent_msg","agent":"Grok","text":"On it. I'll use stone.","actions":[{"type":"place_block","x":10,"y":0,"z":5,"block":"stone"}]}
{"ts":"2026-03-21T14:30:05Z","type":"player_msg","agent":"Grok","text":"Make it taller"}
```

**Agent-World Interaction Log** (`cache/worlds/{id}/ai/world_log.jsonl`):
```jsonl
{"ts":"2026-03-21T14:30:01Z","agent":"Grok","action":"place_block","x":10,"y":0,"z":5,"block":"stone"}
{"ts":"2026-03-21T14:30:02Z","agent":"Grok","action":"move_to","x":11,"y":0,"z":5}
{"ts":"2026-03-21T14:30:03Z","agent":"Grok","action":"place_block","x":11,"y":0,"z":5,"block":"stone"}
{"ts":"2026-03-21T14:30:10Z","agent":"Grok","action":"sleep","reason":"player_out_of_range"}
```

**What we can analyze from logs:**
- Conversation quality — are responses coherent? Does the agent understand requests?
- Action reliability — how often does the agent produce valid vs invalid actions?
- Build patterns — what does the agent actually construct? Does creativity slider change output?
- Performance — response latency, tokens/sec, action throughput
- Memory effectiveness — does the agent remember past conversations and projects?

#### 3d: Memory + Project Persistence
- SQLite-backed conversation memory (`ai_memory.c`)
- Periodic conversation summarization (using the same LLM)
- Project tracking — agent remembers what it was building, can resume after sleep/wake
- Per-world, per-agent memory isolation (`cache/worlds/{id}/ai/{agent_name}/`)

**Deliverable:** Agent walks around the world, responds to commands, builds structures, and all interactions are logged to JSONL files for later analysis.

---

### Future Phases (Not Scoped Yet)

- **Wanderer + Villager classes** — autonomous roaming, trade, self-directed settlement building
- **Inter-agent communication** — villagers coordinate builds, wanderers report discoveries
- **Voice input** — whisper.cpp integration, push-to-talk
- **Web/Mobile AI** — cloud proxy approach TBD, scripted fallback for now
- **Background simulation** — agents keep working when player is far away

---

## Key Design Decisions

### Why Sidecar Server (not Embedded Library)?

**Start with sidecar (`llama-server` as child process):**
- Zero build complexity (no linking libllama into CMake)
- Easy model/prompt iteration without recompiling the game
- Clean process isolation (LLM crash doesn't crash the game)
- You already have WinHTTP patterns from `lobby.c`

**Migrate to embedded later if needed:**
- Single executable distribution
- Lower first-token latency (~50ms REST overhead eliminated)
- Full control over threading and memory
- llama.cpp API is unstable — wait for it to stabilize

### Why Not Cloud AI?

- Latency: local inference at 80+ tok/s beats any cloud roundtrip
- Cost: no API fees, runs on player's hardware
- Offline: works without internet
- Privacy: conversations stay local
- Fits the project philosophy: self-contained, no external dependencies

### Why Grammar Constraints Over Fine-Tuning?

- Works on any model without training
- 100% output validity guaranteed at the token level
- Swap models freely without retraining
- Fine-tuning is expensive and locks you to one model family
