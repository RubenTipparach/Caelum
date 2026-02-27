import { getRedis } from "../../lib/redis.js";

const CORS_HEADERS = {
  "Access-Control-Allow-Origin": "*",
  "Access-Control-Allow-Methods": "POST, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

export default async function handler(req, res) {
  if (req.method === "OPTIONS") {
    return res.status(200).json({});
  }
  if (req.method !== "POST") {
    return res.status(405).json({ error: "Method not allowed" });
  }

  Object.entries(CORS_HEADERS).forEach(([k, v]) => res.setHeader(k, v));

  const { lobby_id, player_name } = req.body || {};
  if (!lobby_id || !player_name) {
    return res.status(400).json({ error: "lobby_id and player_name are required" });
  }

  const redis = getRedis();
  const lobby = await redis.hgetall(`lobby:${lobby_id}`);

  if (!lobby || !lobby.host_name) {
    return res.status(404).json({ error: "Lobby not found" });
  }

  if (Number(lobby.player_count) >= Number(lobby.max_players)) {
    return res.status(409).json({ error: "Lobby is full" });
  }

  // Assign client ID and increment player count
  const client_id = Number(lobby.next_client_id);
  await redis.hset(`lobby:${lobby_id}`, {
    player_count: Number(lobby.player_count) + 1,
    next_client_id: client_id + 1,
  });

  return res.status(200).json({ lobby_id, client_id, player_name });
}
