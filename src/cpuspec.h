#ifndef CPUSPEC_H
#define CPUSPEC_H

#if defined (_TOM_AVX) || defined (_TOM_SSE3)
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>

#ifdef __AVX__
#define BYTE_ALIGNMENT 32
#else
#define BYTE_ALIGNMENT 16
#endif

#ifdef _MSC_VER
#define ALIGNED(X) __declspec(align(BYTE_ALIGNMENT)) X
#else
#define ALIGNED(X) X __attribute__ ((aligned(BYTE_ALIGNMENT)))
#endif
#endif

#endif
