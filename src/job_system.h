#ifndef JOB_SYSTEM_H
#define JOB_SYSTEM_H

#include <stdbool.h>

// Job function signature: takes opaque user data
typedef void (*JobFunc)(void* data);

// Opaque job system handle
typedef struct JobSystem JobSystem;

// Create a job system with the specified number of worker threads.
// Workers start immediately and wait for jobs.
JobSystem* job_system_create(int num_workers);

// Submit a job for execution on a worker thread.
// data is passed to func when it runs. Caller owns data lifetime.
// WARNING: Blocks if queue is full — avoid calling from the main/render thread.
void job_system_submit(JobSystem* sys, JobFunc func, void* data);

// Non-blocking submit: returns true if job was enqueued, false if queue is full.
// Safe to call from the main thread (never blocks).
bool job_system_try_submit(JobSystem* sys, JobFunc func, void* data);

// Get the number of pending (unstarted) jobs in the queue (thread-safe).
int job_system_pending(JobSystem* sys);

// Destroy the job system. Waits for all pending jobs to finish.
void job_system_destroy(JobSystem* sys);

#endif
