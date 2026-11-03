#define _GNU_SOURCE
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#if defined(__i386__) || defined(__x86_64__)
#include <xmmintrin.h>
#define DSPEC_PROFILE_X86 1
#else
#define DSPEC_PROFILE_X86 0
#endif

__attribute__((constructor)) static void dspec_profile_enable_ftz_daz(void) {
  const char *request = getenv("DSPEC_PROFILE_FTZ");
  if (request == NULL || strcmp(request, "1") != 0) return;

#if DSPEC_PROFILE_X86
#if defined(__GNUC__) || defined(__clang__)
  __builtin_cpu_init();
  if (!__builtin_cpu_supports("sse") || !__builtin_cpu_supports("sse2")) {
    (void)fprintf(stderr,
        "DSPEC_PROFILE_FTZ=1 requested, but this x86 CPU lacks required SSE/SSE2 support; aborting diagnostic run.\n");
    (void)fflush(stderr);
    _exit(78);
  }
#endif
  const unsigned int old_mxcsr = _mm_getcsr();
  const unsigned int new_mxcsr = old_mxcsr | (1u << 15) | (1u << 6);
  _mm_setcsr(new_mxcsr);
  (void)fprintf(stderr,
      "DSPEC_PROFILE_FTZ enabled: MXCSR 0x%08x -> 0x%08x (FTZ=1, DAZ=1); subnormal handling is changed.\n",
      old_mxcsr, _mm_getcsr());
  (void)fflush(stderr);
#else
  (void)fprintf(stderr,
      "DSPEC_PROFILE_FTZ=1 requested, but this preload diagnostic supports x86 only; aborting.\n");
  (void)fflush(stderr);
  _exit(78);
#endif
}
