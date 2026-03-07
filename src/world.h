#ifndef WORLD_H
#define WORLD_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

// Forward declarations
struct Camera;

#define WORLD_ID_LEN 12      // hex string, e.g. "a1b2c3d4e5f6"
#define MAX_WORLDS   32

typedef struct WorldMeta {
    char     id[WORLD_ID_LEN + 1];
    char     name[64];
    uint64_t created_at;
    uint64_t last_played;
} WorldMeta;

typedef struct WorldList {
    WorldMeta worlds[MAX_WORLDS];
    int       count;
    int       active_index;   // last played, -1 if none
} WorldList;

// ---- World management ----

void world_generate_id(char* out);
bool world_list_load(WorldList* list);
bool world_list_save(const WorldList* list);

// Create a new world (generates ID, names it "World N", creates dirs). Returns index.
int  world_create_new(WorldList* list);

// Delete a world by index (removes files + shifts list). Returns true on success.
bool world_delete(WorldList* list, int index);

// Path helpers — write into caller-provided buffers
void world_get_dir(const WorldMeta* w, char* out, int size);
void world_get_edits_dir(const WorldMeta* w, char* out, int size);
void world_get_player_path(const WorldMeta* w, char* out, int size);

// Migrate legacy save (cache/player.dat + cache/edits/) to new world system
void world_migrate_legacy(WorldList* list);

// ---- World data pack/unpack (multiplayer) ----

#define WORLD_BLOB_MAGIC 0x574F5244  // "WORD"

typedef struct WorldBlobHeader {
    uint32_t magic;
    uint32_t version;
    int32_t  seed;
    double   host_pos[3];
    float    host_yaw, host_pitch;
    int32_t  file_count;
} WorldBlobHeader;

typedef struct EditBlobEntry {
    char     filename[64];
    uint32_t offset;
    uint32_t size;
} EditBlobEntry;

// Pack all edit files + host position into a binary blob. Caller must free().
uint8_t* world_pack(const WorldMeta* world, const struct Camera* cam, size_t* out_size);

// Unpack a blob: create new world dir, write edit files, create player.dat from host pos.
bool world_unpack(const uint8_t* blob, size_t size, WorldList* list, int* out_index);

#endif
