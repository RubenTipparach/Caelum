import { getRedis } from "../../lib/redis.js";

const CORS_HEADERS = {
  "Access-Control-Allow-Origin": "*",
  "Access-Control-Allow-Methods": "POST, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

function generateLobbyId() {
  const chars = "ABCDEFGHJKLMNPQRSTUVWXYZ23456789";
  let id = "";
  for (let i = 0; i < 6; i++) {
    id += chars[Math.floor(Math.random() * chars.length)];
  }
  return id;
}

export default async function handler(req, res) {
  if (req.method === "OPTIONS") {
    return res.status(200).json({});
  }
  if (req.method !== "POST") {
    return res.status(405).json({ error: "Method not allowed" });
  }

  Object.entries(CORS_HEADERS).forEach(([k, v]) => res.setHeader(k, v));

  const { host_name, max_players = 8 } = req.body || {};
  if (!host_name) {
    return res.status(400).json({ error: "host_name is required" });
  }

  const redis = getRedis();
  const lobby_id = generateLobbyId();

  const lobby = {
    host_name,
    player_count: 1,
    max_players: Math.min(Math.max(max_players, 2), 16),
    created_at: Date.now(),
    status: "waiting",
    next_client_id: 1,
  };

  // Store lobby with 30-minute TTL
  await redis.hset(`lobby:${lobby_id}`, lobby);
  await redis.expire(`lobby:${lobby_id}`, 1800);
  await redis.sadd("active_lobbies", lobby_id);

  return res.status(200).json({ lobby_id, host_name });
}
