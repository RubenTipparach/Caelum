#ifdef _WIN32
#define _CRT_RAND_S
#endif

#include "world.h"
#include "camera.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <direct.h>
#pragma warning(disable: 4996)
#else
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>
#endif

// ---- Platform helpers ----

static void ensure_dir(const char* path) {
#ifdef _WIN32
    CreateDirectoryA(path, NULL);
#else
    mkdir(path, 0755);
#endif
}

static void ensure_dir_recursive(const char* path) {
    char tmp[512];
    snprintf(tmp, sizeof(tmp), "%s", path);
    for (char* p = tmp + 1; *p; p++) {
        if (*p == '/' || *p == '\\') {
            *p = '\0';
            ensure_dir(tmp);
            *p = '/';
        }
    }
    ensure_dir(tmp);
}

static bool file_exists(const char* path) {
#ifdef _WIN32
    DWORD attr = GetFileAttributesA(path);
    return attr != INVALID_FILE_ATTRIBUTES;
#else
    return access(path, F_OK) == 0;
#endif
}

static bool dir_exists(const char* path) {
#ifdef _WIN32
    DWORD attr = GetFileAttributesA(path);
    return attr != INVALID_FILE_ATTRIBUTES && (attr & FILE_ATTRIBUTE_DIRECTORY);
#else
    struct stat st;
    return stat(path, &st) == 0 && S_ISDIR(st.st_mode);
#endif
}

// Recursively delete a directory and all its contents
static void remove_dir_recursive(const char* path) {
#ifdef _WIN32
    WIN32_FIND_DATAA fd;
    char pattern[512];
    snprintf(pattern, sizeof(pattern), "%s\\*", path);
    HANDLE h = FindFirstFileA(pattern, &fd);
    if (h != INVALID_HANDLE_VALUE) {
        do {
            if (fd.cFileName[0] == '.' && (fd.cFileName[1] == '\0' ||
                (fd.cFileName[1] == '.' && fd.cFileName[2] == '\0')))
                continue;
            char child[512];
            snprintf(child, sizeof(child), "%s\\%s", path, fd.cFileName);
            if (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
                remove_dir_recursive(child);
            else
                DeleteFileA(child);
        } while (FindNextFileA(h, &fd));
        FindClose(h);
    }
    RemoveDirectoryA(path);
#else
    DIR* d = opendir(path);
    if (!d) return;
    struct dirent* ent;
    while ((ent = readdir(d)) != NULL) {
        if (ent->d_name[0] == '.' && (ent->d_name[1] == '\0' ||
            (ent->d_name[1] == '.' && ent->d_name[2] == '\0')))
            continue;
        char child[512];
        snprintf(child, sizeof(child), "%s/%s", path, ent->d_name);
        struct stat st;
        if (stat(child, &st) == 0 && S_ISDIR(st.st_mode))
            remove_dir_recursive(child);
        else
            remove(child);
    }
    closedir(d);
    rmdir(path);
#endif
}

// Copy a file from src to dst
static bool copy_file(const char* src, const char* dst) {
#ifdef _WIN32
    return CopyFileA(src, dst, FALSE) != 0;
#else
    FILE* in = fopen(src, "rb");
    if (!in) return false;
    FILE* out = fopen(dst, "wb");
    if (!out) { fclose(in); return false; }
    char buf[4096];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), in)) > 0)
        fwrite(buf, 1, n, out);
    fclose(in);
    fclose(out);
    return true;
#endif
}

// Move (rename) a file
static bool move_file(const char* src, const char* dst) {
#ifdef _WIN32
    return MoveFileA(src, dst) != 0;
#else
    return rename(src, dst) == 0;
#endif
}

// ---- World ID generation ----

void world_generate_id(char* out) {
    static const char hex[] = "0123456789abcdef";

#ifdef _WIN32
    for (int i = 0; i < WORLD_ID_LEN; i++) {
        unsigned int val;
        rand_s(&val);
        out[i] = hex[val % 16];
    }
#else
    FILE* f = fopen("/dev/urandom", "rb");
    if (f) {
        unsigned char bytes[WORLD_ID_LEN];
        fread(bytes, 1, WORLD_ID_LEN, f);
        fclose(f);
        for (int i = 0; i < WORLD_ID_LEN; i++)
            out[i] = hex[bytes[i] % 16];
    } else {
        srand((unsigned)time(NULL));
        for (int i = 0; i < WORLD_ID_LEN; i++)
            out[i] = hex[rand() % 16];
    }
#endif
    out[WORLD_ID_LEN] = '\0';
}

// ---- World list persistence ----

#define WORLDS_DAT_MAGIC 0x574C5354  // "WLST"
#define WORLDS_DAT_VERSION 1

typedef struct WorldsFileHeader {
    uint32_t magic;
    uint32_t version;
    int32_t  count;
    int32_t  active_index;
} WorldsFileHeader;

bool world_list_load(WorldList* list) {
    memset(list, 0, sizeof(WorldList));
    list->active_index = -1;

    FILE* f = fopen("cache/worlds.dat", "rb");
    if (!f) return false;

    WorldsFileHeader hdr;
    if (fread(&hdr, sizeof(hdr), 1, f) != 1 || hdr.magic != WORLDS_DAT_MAGIC) {
        fclose(f);
        return false;
    }

    list->count = hdr.count;
    if (list->count > MAX_WORLDS) list->count = MAX_WORLDS;
    list->active_index = hdr.active_index;

    for (int i = 0; i < list->count; i++) {
        if (fread(&list->worlds[i], sizeof(WorldMeta), 1, f) != 1) {
            list->count = i;
            break;
        }
    }

    fclose(f);
    printf("[WORLD] Loaded %d worlds from cache/worlds.dat\n", list->count);
    fflush(stdout);
    return true;
}

bool world_list_save(const WorldList* list) {
    ensure_dir("cache");

    FILE* f = fopen("cache/worlds.dat", "wb");
    if (!f) {
        printf("[WORLD] ERROR: cannot write cache/worlds.dat\n");
        fflush(stdout);
        return false;
    }

    WorldsFileHeader hdr = {
        .magic = WORLDS_DAT_MAGIC,
        .version = WORLDS_DAT_VERSION,
        .count = list->count,
        .active_index = list->active_index,
    };
    fwrite(&hdr, sizeof(hdr), 1, f);

    for (int i = 0; i < list->count; i++) {
        fwrite(&list->worlds[i], sizeof(WorldMeta), 1, f);
    }

    fclose(f);
    return true;
}

// ---- Path helpers ----

void world_get_dir(const WorldMeta* w, char* out, int size) {
    snprintf(out, size, "cache/worlds/%s", w->id);
}

void world_get_edits_dir(const WorldMeta* w, char* out, int size) {
    snprintf(out, size, "cache/worlds/%s/edits", w->id);
}

void world_get_player_path(const WorldMeta* w, char* out, int size) {
    snprintf(out, size, "cache/worlds/%s/player.dat", w->id);
}

// ---- Create new world ----

int world_create_new(WorldList* list) {
    if (list->count >= MAX_WORLDS) {
        printf("[WORLD] ERROR: max worlds (%d) reached\n", MAX_WORLDS);
        fflush(stdout);
        return -1;
    }

    int idx = list->count;
    WorldMeta* w = &list->worlds[idx];
    memset(w, 0, sizeof(WorldMeta));

    world_generate_id(w->id);
    snprintf(w->name, sizeof(w->name), "World %d", idx + 1);
    w->created_at = (uint64_t)time(NULL);
    w->last_played = w->created_at;

    // Create directories
    char dir[256];
    world_get_dir(w, dir, sizeof(dir));
    ensure_dir_recursive(dir);

    char edits_dir[256];
    world_get_edits_dir(w, edits_dir, sizeof(edits_dir));
    ensure_dir_recursive(edits_dir);

    // Create seed subdirectory inside edits
    char seed_dir[256];
    snprintf(seed_dir, sizeof(seed_dir), "%s/42", edits_dir);
    ensure_dir(seed_dir);

    list->count++;
    list->active_index = idx;

    printf("[WORLD] Created world '%s' (id=%s) at %s\n", w->name, w->id, dir);
    fflush(stdout);

    world_list_save(list);
    return idx;
}

// ---- Delete world ----

bool world_delete(WorldList* list, int index) {
    if (index < 0 || index >= list->count) return false;

    WorldMeta* w = &list->worlds[index];
    printf("[WORLD] Deleting world '%s' (id=%s)\n", w->name, w->id);
    fflush(stdout);

    // Delete world directory
    char dir[256];
    world_get_dir(w, dir, sizeof(dir));
    if (dir_exists(dir))
        remove_dir_recursive(dir);

    // Shift remaining worlds down
    for (int i = index; i < list->count - 1; i++)
        list->worlds[i] = list->worlds[i + 1];
    list->count--;

    // Fix active_index
    if (list->active_index == index)
        list->active_index = -1;
    else if (list->active_index > index)
        list->active_index--;

    world_list_save(list);
    return true;
}

// ---- Legacy migration ----

void world_migrate_legacy(WorldList* list) {
    // Only migrate if there's no worlds.dat and legacy files exist
    if (file_exists("cache/worlds.dat")) return;
    if (!file_exists("cache/player.dat") && !dir_exists("cache/edits")) return;

    printf("[WORLD] Migrating legacy save data...\n");
    fflush(stdout);

    // Create a new world entry
    int idx = world_create_new(list);
    if (idx < 0) return;

    WorldMeta* w = &list->worlds[idx];
    snprintf(w->name, sizeof(w->name), "Legacy World");

    // Move player.dat
    if (file_exists("cache/player.dat")) {
        char player_path[256];
        world_get_player_path(w, player_path, sizeof(player_path));
        if (move_file("cache/player.dat", player_path)) {
            printf("[WORLD] Migrated player.dat -> %s\n", player_path);
        } else {
            // Fallback: copy
            copy_file("cache/player.dat", player_path);
            printf("[WORLD] Copied player.dat -> %s\n", player_path);
        }
    }

    // Move edits directory contents
    if (dir_exists("cache/edits")) {
        char edits_dir[256];
        world_get_edits_dir(w, edits_dir, sizeof(edits_dir));

        // Move cache/edits/42/*.bin -> cache/worlds/{id}/edits/42/*.bin
        char src_seed_dir[256];
        snprintf(src_seed_dir, sizeof(src_seed_dir), "cache/edits/42");

        if (dir_exists(src_seed_dir)) {
            char dst_seed_dir[256];
            snprintf(dst_seed_dir, sizeof(dst_seed_dir), "%s/42", edits_dir);
            ensure_dir_recursive(dst_seed_dir);

#ifdef _WIN32
            WIN32_FIND_DATAA fd;
            char pattern[256];
            snprintf(pattern, sizeof(pattern), "%s/*.bin", src_seed_dir);
            HANDLE hFind = FindFirstFileA(pattern, &fd);
            if (hFind != INVALID_HANDLE_VALUE) {
                do {
                    char src[512], dst[512];
                    snprintf(src, sizeof(src), "%s/%s", src_seed_dir, fd.cFileName);
                    snprintf(dst, sizeof(dst), "%s/%s", dst_seed_dir, fd.cFileName);
                    move_file(src, dst);
                } while (FindNextFileA(hFind, &fd));
                FindClose(hFind);
            }
#else
            DIR* d = opendir(src_seed_dir);
            if (d) {
                struct dirent* ent;
                while ((ent = readdir(d)) != NULL) {
                    if (strstr(ent->d_name, ".bin")) {
                        char src[512], dst[512];
                        snprintf(src, sizeof(src), "%s/%s", src_seed_dir, ent->d_name);
                        snprintf(dst, sizeof(dst), "%s/%s", dst_seed_dir, ent->d_name);
                        move_file(src, dst);
                    }
                }
                closedir(d);
            }
#endif
            printf("[WORLD] Migrated edit sectors -> %s\n", dst_seed_dir);
        }
    }

    fflush(stdout);
    world_list_save(list);
}

// ---- World data pack/unpack (for multiplayer) ----

uint8_t* world_pack(const WorldMeta* world, const struct Camera* cam, size_t* out_size) {
    *out_size = 0;

    // Enumerate edit files
    char edits_dir[256];
    world_get_edits_dir(world, edits_dir, sizeof(edits_dir));
    char seed_dir[256];
    snprintf(seed_dir, sizeof(seed_dir), "%s/42", edits_dir);

    // Collect file names and data
    typedef struct {
        char name[64];
        uint8_t* data;
        uint32_t size;
    } FileEntry;

    FileEntry files[256];
    int file_count = 0;

#ifdef _WIN32
    WIN32_FIND_DATAA fd;
    char pattern[256];
    snprintf(pattern, sizeof(pattern), "%s/*.bin", seed_dir);
    HANDLE hFind = FindFirstFileA(pattern, &fd);
    if (hFind != INVALID_HANDLE_VALUE) {
        do {
            if (file_count >= 256) break;
            char path[512];
            snprintf(path, sizeof(path), "%s/%s", seed_dir, fd.cFileName);
            FILE* f = fopen(path, "rb");
            if (f) {
                fseek(f, 0, SEEK_END);
                long sz = ftell(f);
                fseek(f, 0, SEEK_SET);
                if (sz > 0 && sz < 1024 * 1024) {
                    files[file_count].data = (uint8_t*)malloc(sz);
                    files[file_count].size = (uint32_t)sz;
                    fread(files[file_count].data, 1, sz, f);
                    snprintf(files[file_count].name, 64, "%s", fd.cFileName);
                    file_count++;
                }
                fclose(f);
            }
        } while (FindNextFileA(hFind, &fd));
        FindClose(hFind);
    }
#else
    DIR* d = opendir(seed_dir);
    if (d) {
        struct dirent* ent;
        while ((ent = readdir(d)) != NULL) {
            if (!strstr(ent->d_name, ".bin")) continue;
            if (file_count >= 256) break;
            char path[512];
            snprintf(path, sizeof(path), "%s/%s", seed_dir, ent->d_name);
            FILE* f = fopen(path, "rb");
            if (f) {
                fseek(f, 0, SEEK_END);
                long sz = ftell(f);
                fseek(f, 0, SEEK_SET);
                if (sz > 0 && sz < 1024 * 1024) {
                    files[file_count].data = (uint8_t*)malloc(sz);
                    files[file_count].size = (uint32_t)sz;
                    fread(files[file_count].data, 1, sz, f);
                    snprintf(files[file_count].name, 64, "%s", ent->d_name);
                    file_count++;
                }
                fclose(f);
            }
        }
        closedir(d);
    }
#endif

    // Calculate total blob size
    size_t total = sizeof(WorldBlobHeader)
                 + (size_t)file_count * sizeof(EditBlobEntry);
    for (int i = 0; i < file_count; i++)
        total += files[i].size;

    uint8_t* blob = (uint8_t*)malloc(total);
    if (!blob) {
        for (int i = 0; i < file_count; i++) free(files[i].data);
        return NULL;
    }

    // Write header
    WorldBlobHeader* hdr = (WorldBlobHeader*)blob;
    hdr->magic = WORLD_BLOB_MAGIC;
    hdr->version = 1;
    hdr->seed = 42;
    hdr->host_pos[0] = cam->pos_d[0];
    hdr->host_pos[1] = cam->pos_d[1];
    hdr->host_pos[2] = cam->pos_d[2];
    hdr->host_yaw = cam->yaw;
    hdr->host_pitch = cam->pitch;
    hdr->file_count = file_count;

    // Write entries + file data
    EditBlobEntry* entries = (EditBlobEntry*)(blob + sizeof(WorldBlobHeader));
    uint32_t data_offset = (uint32_t)(sizeof(WorldBlobHeader) + file_count * sizeof(EditBlobEntry));

    for (int i = 0; i < file_count; i++) {
        snprintf(entries[i].filename, 64, "%s", files[i].name);
        entries[i].offset = data_offset;
        entries[i].size = files[i].size;
        memcpy(blob + data_offset, files[i].data, files[i].size);
        data_offset += files[i].size;
        free(files[i].data);
    }

    *out_size = total;
    printf("[WORLD] Packed %d edit files, %zu bytes total\n", file_count, total);
    fflush(stdout);
    return blob;
}

bool world_unpack(const uint8_t* blob, size_t size, WorldList* list, int* out_index) {
    if (size < sizeof(WorldBlobHeader)) return false;

    const WorldBlobHeader* hdr = (const WorldBlobHeader*)blob;
    if (hdr->magic != WORLD_BLOB_MAGIC || hdr->version != 1) return false;

    size_t expected_min = sizeof(WorldBlobHeader) + (size_t)hdr->file_count * sizeof(EditBlobEntry);
    if (size < expected_min) return false;

    // Create new world for the client
    int idx = world_create_new(list);
    if (idx < 0) return false;

    WorldMeta* w = &list->worlds[idx];
    snprintf(w->name, sizeof(w->name), "MP World %d", idx + 1);

    // Write edit files
    char edits_dir[256];
    world_get_edits_dir(w, edits_dir, sizeof(edits_dir));
    char seed_dir[256];
    snprintf(seed_dir, sizeof(seed_dir), "%s/%d", edits_dir, hdr->seed);
    ensure_dir_recursive(seed_dir);

    const EditBlobEntry* entries = (const EditBlobEntry*)(blob + sizeof(WorldBlobHeader));
    for (int i = 0; i < hdr->file_count; i++) {
        if (entries[i].offset + entries[i].size > size) continue;
        char path[512];
        snprintf(path, sizeof(path), "%s/%s", seed_dir, entries[i].filename);
        FILE* f = fopen(path, "wb");
        if (f) {
            fwrite(blob + entries[i].offset, 1, entries[i].size, f);
            fclose(f);
        }
    }

    // Write player.dat from host position
    char player_path[256];
    world_get_player_path(w, player_path, sizeof(player_path));
    FILE* f = fopen(player_path, "wb");
    if (f) {
        uint32_t magic = 0x504C5952;  // PLYR
        fwrite(&magic, sizeof(magic), 1, f);
        fwrite(hdr->host_pos, sizeof(double), 3, f);
        fwrite(&hdr->host_yaw, sizeof(float), 1, f);
        fwrite(&hdr->host_pitch, sizeof(float), 1, f);
        uint8_t slot = 0;
        fwrite(&slot, 1, 1, f);
        uint8_t pad[3] = {0};
        fwrite(pad, 1, 3, f);
        fclose(f);
    }

    world_list_save(list);
    if (out_index) *out_index = idx;

    printf("[WORLD] Unpacked world '%s' with %d edit files\n", w->name, hdr->file_count);
    fflush(stdout);
    return true;
}
