import { getRedis } from "../../lib/redis.js";

const CORS_HEADERS = {
  "Access-Control-Allow-Origin": "*",
  "Access-Control-Allow-Methods": "GET, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

export default async function handler(req, res) {
  if (req.method === "OPTIONS") {
    return res.status(200).json({});
  }
  if (req.method !== "GET") {
    return res.status(405).json({ error: "Method not allowed" });
  }

  Object.entries(CORS_HEADERS).forEach(([k, v]) => res.setHeader(k, v));

  const redis = getRedis();
  const lobby_ids = await redis.smembers("active_lobbies");

  const lobbies = [];
  for (const id of lobby_ids) {
    const lobby = await redis.hgetall(`lobby:${id}`);
    if (!lobby || !lobby.host_name) {
      // Expired lobby — clean up from the set
      await redis.srem("active_lobbies", id);
      continue;
    }
    if (Number(lobby.player_count) >= Number(lobby.max_players)) {
      continue; // Full lobby
    }
    lobbies.push({
      lobby_id: id,
      host_name: lobby.host_name,
      player_count: Number(lobby.player_count),
      max_players: Number(lobby.max_players),
      created_at: Number(lobby.created_at),
    });
  }

  return res.status(200).json({ lobbies });
}
