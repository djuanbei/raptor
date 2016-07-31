#ifndef __SYSTEM_H
#define __SYSTEM_H

#if defined(__linux__)
#include <fpu_control.h>
#endif

#include "IntTypes.h"

//-------------------------------------------------------------------------------------------------
static inline double systemTime( void );// SYSTEM time in seconds.
static inline double cpuTime(void);  // CPU-time in seconds.

extern double memUsed();  // Memory in mega bytes (returns 0 for unsupported
                          // architectures).
extern double memUsedPeak(bool strictlyPeak = false);  // Peak-memory in mega
                                                       // bytes (returns 0 for
                                                       // unsupported
                                                       // architectures).

extern void setX86FPUPrecision();  // Make sure double's are represented with
                                   // the same precision
                                   // in memory and registers.

extern void limitMemory(
    uint64_t max_mem_mb);  // Set a limit on total memory usage. The exact
                           // semantics varies depending on architecture.

extern void limitTime(
    uint32_t max_cpu_time);  // Set a limit on maximum CPU time. The exact
                             // semantics varies depending on architecture.

extern void sigTerm(
    void handler(int));  // Set up handling of available termination signals.

//-------------------------------------------------------------------------------------------------
// Implementation of inline functions:

#if defined(_MSC_VER) || defined(__MINGW32__)
#include <time.h>

static inline double cpuTime(void) { return (double)clock() / CLOCKS_PER_SEC; }

#else
#include<time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

static inline double cpuTime(void) {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

static inline double systemTime( void ){
  struct timespec start;
  clock_gettime( CLOCK_MONOTONIC, &start);

  return start.tv_sec +start.tv_nsec/1000000000.0;
      
}

#endif

#endif
