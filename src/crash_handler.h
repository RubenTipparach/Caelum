#ifndef CRASH_HANDLER_H
#define CRASH_HANDLER_H

// Install a crash handler that writes diagnostic info to crash.log
// Call this as early as possible in main/init.
void crash_handler_install(void);

// Sokol-compatible logger that writes errors to crash.log (and also to console).
// Pass as .logger.func in sg_setup, sdtx_setup, etc.
void crash_log_func(const char* tag, unsigned int log_level, unsigned int log_item,
                    const char* message, unsigned int line_nr,
                    const char* filename, void* user_data);

#endif
