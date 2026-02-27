import { getRedis } from "../../lib/redis.js";

const CORS_HEADERS = {
  "Access-Control-Allow-Origin": "*",
  "Access-Control-Allow-Methods": "GET, POST, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

export default async function handler(req, res) {
  Object.entries(CORS_HEADERS).forEach(([k, v]) => res.setHeader(k, v));

  if (req.method === "OPTIONS") {
    return res.status(200).json({});
  }

  const redis = getRedis();

  if (req.method === "POST") {
    // Broadcast a block edit to all players in the lobby
    const { lobby_id, player_id, cell_index, action, color } = req.body || {};
    if (!lobby_id || player_id === undefined || cell_index === undefined || !action) {
      return res.status(400).json({
        error: "lobby_id, player_id, cell_index, and action are required",
      });
    }

    const edit = {
      player_id,
      cell_index,
      action, // "break" or "place"
      color: color || null,
      timestamp: Date.now(),
    };

    // Push edit to lobby's edit stream (all players poll this)
    const key = `edits:${lobby_id}`;
    await redis.lpush(key, JSON.stringify(edit));
    await redis.expire(key, 1800); // 30 min TTL

    // Keep only last 1000 edits in the stream
    await redis.ltrim(key, 0, 999);

    return res.status(200).json({ ok: true });
  }

  if (req.method === "GET") {
    // Poll for new edits since a given timestamp
    const { lobby_id, since } = req.query || {};
    if (!lobby_id) {
      return res.status(400).json({ error: "lobby_id is required" });
    }

    const key = `edits:${lobby_id}`;
    const sinceTs = since ? Number(since) : 0;

    // Get all edits (newest first)
    const raw = await redis.lrange(key, 0, -1);
    const edits = [];
    for (const item of raw) {
      try {
        const edit = typeof item === "string" ? JSON.parse(item) : item;
        if (edit.timestamp > sinceTs) {
          edits.push(edit);
        }
      } catch {
        // skip malformed
      }
    }

    // Return in chronological order (oldest first)
    edits.reverse();

    return res.status(200).json({ edits });
  }

  return res.status(405).json({ error: "Method not allowed" });
}
