#include <signal.h>
#include <stdio.h>

#include "System.h"


#if defined(__linux__)

#include <stdlib.h>



static inline int memReadStat(int field)
{
    char  name[256];
    pid_t pid = getpid();
    int   value;

    sprintf(name, "/proc/%d/statm", pid);
    FILE* in = fopen(name, "rb");
    if (in == NULL) return 0;

    for (; field >= 0; field--)
        if (fscanf(in, "%d", &value) != 1)
            printf("ERROR! Failed to parse memory statistics from \"/proc\".\n"), exit(1);
    fclose(in);
    return value;
}


static inline int memReadPeak(void)
{
    char  name[256];
    pid_t pid = getpid();

    sprintf(name, "/proc/%d/status", pid);
    FILE* in = fopen(name, "rb");
    if (in == NULL) return 0;

    // Find the correct line, beginning with "VmPeak:":
    int peak_kb = 0;
    while (!feof(in) && fscanf(in, "VmPeak: %d kB", &peak_kb) != 1)
        while (!feof(in) && fgetc(in) != '\n')
            ;
    fclose(in);

    return peak_kb;
}

double memUsed() { return (double)memReadStat(0) * (double)getpagesize() / (1024*1024); }
double memUsedPeak(bool strictlyPeak) { 
    double peak = memReadPeak() / (double)1024;
    return peak == 0 && !strictlyPeak ? memUsed() : peak; }

#elif defined(__FreeBSD__) || defined(__FreeBSD_kernel__) || defined(__gnu_hurd__)

double memUsed() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_maxrss / 1024; }
double memUsedPeak() { return memUsed(); }


#elif defined(__APPLE__)
#include <malloc/malloc.h>

double memUsed() {
    malloc_statistics_t t;
    malloc_zone_statistics(NULL, &t);
    return (double)t.max_size_in_use / (1024*1024); }
double memUsedPeak() { return memUsed(); }

#else
double memUsed()     { return 0; }
double memUsedPeak() { return 0; }
#endif


void setX86FPUPrecision()
{
#if defined(__linux__) && defined(_FPU_EXTENDED) && defined(_FPU_DOUBLE) && defined(_FPU_GETCW)
    // Only correct FPU precision on Linux architectures that needs and supports it:
    fpu_control_t oldcw, newcw;
    _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
    printf("WARNING: for repeatability, setting FPU to use double precision\n");
#endif
}


#if !defined(_MSC_VER) && !defined(__MINGW32__)
void limitMemory(uint64_t max_mem_mb)
{
// FIXME: OpenBSD does not support RLIMIT_AS. Not sure how well RLIMIT_DATA works instead.
#if defined(__OpenBSD__)
#define RLIMIT_AS RLIMIT_DATA
#endif

    // Set limit on virtual memory:
    if (max_mem_mb != 0){
        rlim_t new_mem_lim = (rlim_t)max_mem_mb * 1024*1024;
        rlimit rl;
        getrlimit(RLIMIT_AS, &rl);
        if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
            rl.rlim_cur = new_mem_lim;
            if (setrlimit(RLIMIT_AS, &rl) == -1)
                printf("WARNING! Could not set resource limit: Virtual memory.\n");
        }
    }

#if defined(__OpenBSD__)
#undef RLIMIT_AS
#endif
}
#else
void limitMemory(uint64_t /*max_mem_mb*/)
{
    printf("WARNING! Memory limit not supported on this architecture.\n");
}
#endif


#if !defined(_MSC_VER) && !defined(__MINGW32__)
void limitTime(uint32_t max_cpu_time)
{
    if (max_cpu_time != 0){
        rlimit rl;
        getrlimit(RLIMIT_CPU, &rl);
        if (rl.rlim_max == RLIM_INFINITY || (rlim_t)max_cpu_time < rl.rlim_max){
            rl.rlim_cur = max_cpu_time;
            if (setrlimit(RLIMIT_CPU, &rl) == -1)
                printf("WARNING! Could not set resource limit: CPU-time.\n");
        }
    }
}
#else
void limitTime(uint32_t /*max_cpu_time*/)
{
    printf("WARNING! CPU-time limit not supported on this architecture.\n");
}
#endif


void sigTerm(void handler(int))
{
    signal(SIGINT, handler);
    signal(SIGTERM,handler);
#ifdef SIGXCPU
    signal(SIGXCPU,handler);
#endif
}