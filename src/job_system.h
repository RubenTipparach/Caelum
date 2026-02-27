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
void job_system_submit(JobSystem* sys, JobFunc func, void* data);

// Get the number of pending (unstarted) jobs
int job_system_pending(const JobSystem* sys);

// Destroy the job system. Waits for all pending jobs to finish.
void job_system_destroy(JobSystem* sys);

#endif
