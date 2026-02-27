#include "crash_handler.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <dbghelp.h>
#pragma comment(lib, "dbghelp.lib")

static const char* exception_code_str(DWORD code) {
    switch (code) {
        case EXCEPTION_ACCESS_VIOLATION:    return "ACCESS_VIOLATION";
        case EXCEPTION_STACK_OVERFLOW:      return "STACK_OVERFLOW";
        case EXCEPTION_INT_DIVIDE_BY_ZERO:  return "INT_DIVIDE_BY_ZERO";
        case EXCEPTION_FLT_DIVIDE_BY_ZERO:  return "FLT_DIVIDE_BY_ZERO";
        case EXCEPTION_ILLEGAL_INSTRUCTION: return "ILLEGAL_INSTRUCTION";
        case EXCEPTION_IN_PAGE_ERROR:       return "IN_PAGE_ERROR";
        case EXCEPTION_ARRAY_BOUNDS_EXCEEDED: return "ARRAY_BOUNDS_EXCEEDED";
        default: return "UNKNOWN";
    }
}

static LONG WINAPI crash_exception_filter(EXCEPTION_POINTERS* ep) {
    FILE* f = fopen("crash.log", "a");
    if (!f) f = stderr;

    time_t now = time(NULL);
    fprintf(f, "\n========== CRASH REPORT ==========\n");
    fprintf(f, "Time: %s", ctime(&now));
    fprintf(f, "Exception: 0x%08lX (%s)\n",
            ep->ExceptionRecord->ExceptionCode,
            exception_code_str(ep->ExceptionRecord->ExceptionCode));
    fprintf(f, "Address: 0x%p\n", ep->ExceptionRecord->ExceptionAddress);

    if (ep->ExceptionRecord->ExceptionCode == EXCEPTION_ACCESS_VIOLATION &&
        ep->ExceptionRecord->NumberParameters >= 2) {
        const char* op = ep->ExceptionRecord->ExceptionInformation[0] == 0 ? "read" : "write";
        fprintf(f, "Fault: %s of address 0x%p\n", op,
                (void*)(uintptr_t)ep->ExceptionRecord->ExceptionInformation[1]);
    }

#ifdef _M_X64
    CONTEXT* ctx = ep->ContextRecord;
    fprintf(f, "\nRegisters:\n");
    fprintf(f, "  RIP = 0x%016llX\n", ctx->Rip);
    fprintf(f, "  RSP = 0x%016llX\n", ctx->Rsp);
    fprintf(f, "  RBP = 0x%016llX\n", ctx->Rbp);
    fprintf(f, "  RAX = 0x%016llX\n", ctx->Rax);
    fprintf(f, "  RBX = 0x%016llX\n", ctx->Rbx);
    fprintf(f, "  RCX = 0x%016llX\n", ctx->Rcx);
    fprintf(f, "  RDX = 0x%016llX\n", ctx->Rdx);
#endif

    // Stack trace using StackWalk64
    fprintf(f, "\nStack trace:\n");
    HANDLE process = GetCurrentProcess();
    HANDLE thread = GetCurrentThread();
    SymInitialize(process, NULL, TRUE);

    STACKFRAME64 frame = {0};
#ifdef _M_X64
    DWORD machine_type = IMAGE_FILE_MACHINE_AMD64;
    frame.AddrPC.Offset = ep->ContextRecord->Rip;
    frame.AddrFrame.Offset = ep->ContextRecord->Rbp;
    frame.AddrStack.Offset = ep->ContextRecord->Rsp;
#else
    DWORD machine_type = IMAGE_FILE_MACHINE_I386;
    frame.AddrPC.Offset = ep->ContextRecord->Eip;
    frame.AddrFrame.Offset = ep->ContextRecord->Ebp;
    frame.AddrStack.Offset = ep->ContextRecord->Esp;
#endif
    frame.AddrPC.Mode = AddrModeFlat;
    frame.AddrFrame.Mode = AddrModeFlat;
    frame.AddrStack.Mode = AddrModeFlat;

    CONTEXT ctx_copy = *ep->ContextRecord;
    char sym_buf[sizeof(SYMBOL_INFO) + 256];

    for (int i = 0; i < 32; i++) {
        if (!StackWalk64(machine_type, process, thread, &frame, &ctx_copy,
                         NULL, SymFunctionTableAccess64, SymGetModuleBase64, NULL)) {
            break;
        }
        if (frame.AddrPC.Offset == 0) break;

        SYMBOL_INFO* sym = (SYMBOL_INFO*)sym_buf;
        sym->SizeOfStruct = sizeof(SYMBOL_INFO);
        sym->MaxNameLen = 255;

        DWORD64 displacement = 0;
        if (SymFromAddr(process, frame.AddrPC.Offset, &displacement, sym)) {
            fprintf(f, "  [%2d] %s + 0x%llX (0x%016llX)\n",
                    i, sym->Name, displacement, frame.AddrPC.Offset);
        } else {
            fprintf(f, "  [%2d] 0x%016llX\n", i, frame.AddrPC.Offset);
        }
    }

    SymCleanup(process);
    fprintf(f, "==================================\n");

    if (f != stderr) {
        fflush(f);
        fclose(f);
    }

    // Also dump to stderr so it shows in console
    fprintf(stderr, "\n[CRASH] %s at 0x%p — see crash.log for details\n",
            exception_code_str(ep->ExceptionRecord->ExceptionCode),
            ep->ExceptionRecord->ExceptionAddress);
    fflush(stderr);

    return EXCEPTION_EXECUTE_HANDLER;
}

void crash_handler_install(void) {
    SetUnhandledExceptionFilter(crash_exception_filter);
    // Clear crash.log at start of each run
    FILE* f = fopen("crash.log", "w");
    if (f) fclose(f);
}

#else
// POSIX: use signal handlers
#include <signal.h>
#include <string.h>
#include <execinfo.h>

static void crash_signal_handler(int sig) {
    FILE* f = fopen("crash.log", "a");
    if (!f) f = stderr;

    time_t now = time(NULL);
    fprintf(f, "\n========== CRASH REPORT ==========\n");
    fprintf(f, "Time: %s", ctime(&now));
    fprintf(f, "Signal: %d (%s)\n", sig, strsignal(sig));

    void* bt[32];
    int bt_size = backtrace(bt, 32);
    char** bt_syms = backtrace_symbols(bt, bt_size);
    if (bt_syms) {
        fprintf(f, "\nStack trace:\n");
        for (int i = 0; i < bt_size; i++) {
            fprintf(f, "  [%2d] %s\n", i, bt_syms[i]);
        }
        free(bt_syms);
    }

    fprintf(f, "==================================\n");
    if (f != stderr) {
        fflush(f);
        fclose(f);
    }

    fprintf(stderr, "\n[CRASH] Signal %d (%s) — see crash.log for details\n",
            sig, strsignal(sig));

    _exit(1);
}

void crash_handler_install(void) {
    signal(SIGSEGV, crash_signal_handler);
    signal(SIGABRT, crash_signal_handler);
    signal(SIGFPE, crash_signal_handler);
    signal(SIGILL, crash_signal_handler);
    // Clear crash.log at start of each run
    FILE* f = fopen("crash.log", "w");
    if (f) fclose(f);
}
#endif

// ---- Sokol-compatible logger that writes errors to crash.log ----
// log_level: 0=panic, 1=error, 2=warn, 3=info
void crash_log_func(const char* tag, unsigned int log_level, unsigned int log_item,
                    const char* message, unsigned int line_nr,
                    const char* filename, void* user_data) {
    (void)user_data;
    (void)filename;

    const char* level_str;
    switch (log_level) {
        case 0: level_str = "panic"; break;
        case 1: level_str = "error"; break;
        case 2: level_str = "warn";  break;
        default: level_str = "info"; break;
    }

    // Always print to console
    printf("[%s][%s][id:%u][line:%u] %s\n",
           tag ? tag : "?", level_str, log_item, line_nr,
           message ? message : "");

    // Write errors and panics to crash.log
    if (log_level <= 1) {
        FILE* f = fopen("crash.log", "a");
        if (f) {
            fprintf(f, "[%s][%s][id:%u][line:%u] %s\n",
                    tag ? tag : "?", level_str, log_item, line_nr,
                    message ? message : "");
            fflush(f);
            fclose(f);
        }
    }
}
