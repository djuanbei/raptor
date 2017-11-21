
#ifndef __INTTYPES_H
#define __INTTYPES_H

#ifdef __sun
// Not sure if there are newer versions that support C99 headers. The
// needed features are implemented in the headers below though:

#include <sys/int_fmtio.h>
#include <sys/int_limits.h>
#include <sys/int_types.h>

#else

#include <inttypes.h>
#include <stdint.h>

#endif

#include <climits>

#endif
