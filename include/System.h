#ifndef __SYSTEM_H
#define __SYSTEM_H

#if defined(__linux__)
#include <fpu_control.h>
#endif

#include "myinttype.h"

#define callTime(t, exper)                                \
  auto t0 = std::chrono::system_clock::now();             \
  exper;                                                  \
  auto t1 = std::chrono::system_clock::now();             \
  std::chrono::duration<double, std::milli> ms = t1 - t0; \
  t = ms.count();


#define LOOPEXP(start, end, A, B, op)                         \
  {                                                        \
  int i=start;                                             \
  for(; i+3< end; i+=4){                                   \
    A[i] =A[i] op  B[i];                                        \
    A[i+1] =A[i+1] op B[i+1];                                    \
    A[i+2] =A[i+2] op B[i+2];                                    \
    A[i+3] = A[i+3] op B[i+3];                                    \
  }                                                        \
  for(; i< end; i++) {                                     \
    A[i]= A[i] op B[i];                                        \
  }}                                                       \




//-------------------------------------------------------------------------------------------------
static inline double systemTime(void);  // SYSTEM time in seconds.
static inline double cpuTime(void);     // CPU-time in seconds.

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

#elif (__linux__)
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

static inline double cpuTime(void) {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

static inline double systemTime(void) {
  struct timespec start;
  clock_gettime(CLOCK_MONOTONIC, &start);

  return start.tv_sec + start.tv_nsec / 1000000000.0;
}

#elif (__MACH__)

#include <sys/time.h>
#include <time.h>

#include <mach/clock.h>
#include <mach/mach.h>
static inline double cpuTime(void) { return (double)clock() / CLOCKS_PER_SEC; }

static inline double systemTime(void) {
  struct timespec ts;

  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts.tv_sec = mts.tv_sec;
  ts.tv_nsec = mts.tv_nsec;

  return ts.tv_sec + ts.tv_nsec / 1000000000.0;
}

#endif

#endif
