#include "job_system.h"

#include <stdlib.h>
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>

// ---- Win32 threading (MSVC doesn't fully support C11 threads) ----

#define JOB_QUEUE_SIZE 256

typedef struct Job {
    JobFunc func;
    void* data;
} Job;

struct JobSystem {
    // Ring buffer queue
    Job queue[JOB_QUEUE_SIZE];
    int queue_head;     // Next slot to dequeue from
    int queue_tail;     // Next slot to enqueue to
    int queue_count;    // Number of items in queue

    // Synchronization
    CRITICAL_SECTION mutex;
    CONDITION_VARIABLE cond_not_empty;   // Signal when job available
    CONDITION_VARIABLE cond_not_full;    // Signal when queue has space

    // Worker threads
    HANDLE* threads;
    int num_workers;
    bool shutdown;
};

static DWORD WINAPI worker_thread(LPVOID param) {
    JobSystem* sys = (JobSystem*)param;

    for (;;) {
        Job job;

        EnterCriticalSection(&sys->mutex);

        // Wait for a job or shutdown
        while (sys->queue_count == 0 && !sys->shutdown) {
            SleepConditionVariableCS(&sys->cond_not_empty, &sys->mutex, INFINITE);
        }

        if (sys->shutdown && sys->queue_count == 0) {
            LeaveCriticalSection(&sys->mutex);
            return 0;
        }

        // Dequeue
        job = sys->queue[sys->queue_head];
        sys->queue_head = (sys->queue_head + 1) % JOB_QUEUE_SIZE;
        sys->queue_count--;

        // Signal that there's space in the queue
        WakeConditionVariable(&sys->cond_not_full);

        LeaveCriticalSection(&sys->mutex);

        // Execute the job outside the lock
        if (job.func) {
            job.func(job.data);
        }
    }
}

JobSystem* job_system_create(int num_workers) {
    JobSystem* sys = (JobSystem*)calloc(1, sizeof(JobSystem));
    sys->num_workers = num_workers;
    sys->shutdown = false;
    sys->queue_head = 0;
    sys->queue_tail = 0;
    sys->queue_count = 0;

    InitializeCriticalSection(&sys->mutex);
    InitializeConditionVariable(&sys->cond_not_empty);
    InitializeConditionVariable(&sys->cond_not_full);

    sys->threads = (HANDLE*)calloc(num_workers, sizeof(HANDLE));
    for (int i = 0; i < num_workers; i++) {
        sys->threads[i] = CreateThread(NULL, 0, worker_thread, sys, 0, NULL);
    }

    printf("[JOBS] Created job system with %d workers\n", num_workers);
    fflush(stdout);
    return sys;
}

void job_system_submit(JobSystem* sys, JobFunc func, void* data) {
    EnterCriticalSection(&sys->mutex);

    // Wait if queue is full
    while (sys->queue_count >= JOB_QUEUE_SIZE) {
        SleepConditionVariableCS(&sys->cond_not_full, &sys->mutex, INFINITE);
    }

    // Enqueue
    sys->queue[sys->queue_tail].func = func;
    sys->queue[sys->queue_tail].data = data;
    sys->queue_tail = (sys->queue_tail + 1) % JOB_QUEUE_SIZE;
    sys->queue_count++;

    // Wake a worker
    WakeConditionVariable(&sys->cond_not_empty);

    LeaveCriticalSection(&sys->mutex);
}

int job_system_pending(const JobSystem* sys) {
    return sys->queue_count;
}

void job_system_destroy(JobSystem* sys) {
    if (!sys) return;

    // Signal shutdown
    EnterCriticalSection(&sys->mutex);
    sys->shutdown = true;
    WakeAllConditionVariable(&sys->cond_not_empty);
    LeaveCriticalSection(&sys->mutex);

    // Wait for all workers to finish
    WaitForMultipleObjects(sys->num_workers, sys->threads, TRUE, INFINITE);

    for (int i = 0; i < sys->num_workers; i++) {
        CloseHandle(sys->threads[i]);
    }
    free(sys->threads);

    DeleteCriticalSection(&sys->mutex);
    free(sys);

    printf("[JOBS] Job system destroyed\n");
    fflush(stdout);
}

#else
// ---- POSIX / Emscripten fallback (pthreads) ----

#include <pthread.h>

#define JOB_QUEUE_SIZE 256

typedef struct Job {
    JobFunc func;
    void* data;
} Job;

struct JobSystem {
    Job queue[JOB_QUEUE_SIZE];
    int queue_head;
    int queue_tail;
    int queue_count;

    pthread_mutex_t mutex;
    pthread_cond_t cond_not_empty;
    pthread_cond_t cond_not_full;

    pthread_t* threads;
    int num_workers;
    bool shutdown;
};

static void* worker_thread(void* param) {
    JobSystem* sys = (JobSystem*)param;

    for (;;) {
        Job job;

        pthread_mutex_lock(&sys->mutex);

        while (sys->queue_count == 0 && !sys->shutdown) {
            pthread_cond_wait(&sys->cond_not_empty, &sys->mutex);
        }

        if (sys->shutdown && sys->queue_count == 0) {
            pthread_mutex_unlock(&sys->mutex);
            return NULL;
        }

        job = sys->queue[sys->queue_head];
        sys->queue_head = (sys->queue_head + 1) % JOB_QUEUE_SIZE;
        sys->queue_count--;

        pthread_cond_signal(&sys->cond_not_full);
        pthread_mutex_unlock(&sys->mutex);

        if (job.func) {
            job.func(job.data);
        }
    }
}

JobSystem* job_system_create(int num_workers) {
    JobSystem* sys = (JobSystem*)calloc(1, sizeof(JobSystem));
    sys->num_workers = num_workers;
    sys->shutdown = false;

    pthread_mutex_init(&sys->mutex, NULL);
    pthread_cond_init(&sys->cond_not_empty, NULL);
    pthread_cond_init(&sys->cond_not_full, NULL);

    sys->threads = (pthread_t*)calloc(num_workers, sizeof(pthread_t));
    for (int i = 0; i < num_workers; i++) {
        pthread_create(&sys->threads[i], NULL, worker_thread, sys);
    }

    printf("[JOBS] Created job system with %d workers\n", num_workers);
    fflush(stdout);
    return sys;
}

void job_system_submit(JobSystem* sys, JobFunc func, void* data) {
    pthread_mutex_lock(&sys->mutex);

    while (sys->queue_count >= JOB_QUEUE_SIZE) {
        pthread_cond_wait(&sys->cond_not_full, &sys->mutex);
    }

    sys->queue[sys->queue_tail].func = func;
    sys->queue[sys->queue_tail].data = data;
    sys->queue_tail = (sys->queue_tail + 1) % JOB_QUEUE_SIZE;
    sys->queue_count++;

    pthread_cond_signal(&sys->cond_not_empty);
    pthread_mutex_unlock(&sys->mutex);
}

int job_system_pending(const JobSystem* sys) {
    return sys->queue_count;
}

void job_system_destroy(JobSystem* sys) {
    if (!sys) return;

    pthread_mutex_lock(&sys->mutex);
    sys->shutdown = true;
    pthread_cond_broadcast(&sys->cond_not_empty);
    pthread_mutex_unlock(&sys->mutex);

    for (int i = 0; i < sys->num_workers; i++) {
        pthread_join(sys->threads[i], NULL);
    }
    free(sys->threads);

    pthread_mutex_destroy(&sys->mutex);
    pthread_cond_destroy(&sys->cond_not_empty);
    pthread_cond_destroy(&sys->cond_not_full);
    free(sys);

    printf("[JOBS] Job system destroyed\n");
    fflush(stdout);
}

#endif
