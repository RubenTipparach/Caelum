#ifndef LOBBY_H
#define LOBBY_H

#include <stdbool.h>
#include <stdint.h>

typedef enum {
    LOBBY_IDLE,
    LOBBY_PENDING,
    LOBBY_SUCCESS,
    LOBBY_ERROR,
} LobbyState;

typedef struct LobbyRequest {
    LobbyState state;
    int        http_status;
    char       body[4096];
    int        body_len;
    char       error[128];

    // Parsed response fields (populated on LOBBY_SUCCESS)
    char       lobby_id[8];     // 6-char room code + null
    int        client_id;
} LobbyRequest;

// Initialize lobby system with API base URL (e.g. "https://hex-planets.vercel.app")
void lobby_init(const char* api_base_url);
void lobby_shutdown(void);

// Create a new lobby (non-blocking). Poll with lobby_poll().
void lobby_create(const char* host_name, LobbyRequest* req);

// Join an existing lobby (non-blocking). Poll with lobby_poll().
void lobby_join(const char* lobby_id, const char* player_name, LobbyRequest* req);

// Poll for async HTTP completion. Call each frame while req->state == LOBBY_PENDING.
void lobby_poll(LobbyRequest* req);

#endif
