import { getRedis } from "../lib/redis.js";

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
    // Send a signaling message (SDP offer/answer or ICE candidate)
    const { lobby_id, from_id, to_id, type, payload } = req.body || {};
    if (!lobby_id || from_id === undefined || to_id === undefined || !type || !payload) {
      return res.status(400).json({
        error: "lobby_id, from_id, to_id, type, and payload are required",
      });
    }

    const message = JSON.stringify({ from_id, type, payload, timestamp: Date.now() });
    const key = `signals:${lobby_id}:${to_id}`;

    await redis.lpush(key, message);
    await redis.expire(key, 300); // 5-minute TTL on signal queues

    return res.status(200).json({ ok: true });
  }

  if (req.method === "GET") {
    // Poll for signaling messages
    const { lobby_id, client_id } = req.query || {};
    if (!lobby_id || client_id === undefined) {
      return res.status(400).json({ error: "lobby_id and client_id are required" });
    }

    const key = `signals:${lobby_id}:${client_id}`;

    // Drain all pending messages
    const messages = [];
    while (true) {
      const msg = await redis.rpop(key);
      if (!msg) break;
      try {
        messages.push(typeof msg === "string" ? JSON.parse(msg) : msg);
      } catch {
        messages.push(msg);
      }
    }

    return res.status(200).json({ messages });
  }

  return res.status(405).json({ error: "Method not allowed" });
}
