#!/usr/bin/env bash
#
# Build and run the test suite under a sanitizer.
#
#   scripts/test_sanitized.sh address [extra cmake args...]
#   scripts/test_sanitized.sh undefined
#
# "address" turns on AddressSanitizer and UndefinedBehaviorSanitizer together,
# which is how this library has always been checked: they compose, and the
# combined run costs little more than either alone.
#
# There is deliberately no "thread" mode. ThreadSanitizer cannot validate this
# library at all: the parallelism is OpenMP, GCC's libgomp carries no TSan
# annotations, and every barrier and reduction is therefore reported as a race.
# That is a limitation of the tooling and not a finding, and a job that reports
# hundreds of false positives is worse than no job.
set -euo pipefail

sanitizer="${1:?usage: $0 address|undefined [cmake args...]}"
shift || true

case "${sanitizer}" in
  address)   flags="-fsanitize=address,undefined" ;;
  undefined) flags="-fsanitize=undefined" ;;
  *) echo "unknown sanitizer: ${sanitizer}" >&2; exit 2 ;;
esac

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build="${root}/build-${sanitizer}"

cmake -S "${root}" -B "${build}" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="${flags} -fno-omit-frame-pointer -g" \
  -DGSHTRANS_BUILD_BENCHMARKS=OFF \
  -DGSHTRANS_INSTALL=OFF \
  "$@"

cmake --build "${build}" --parallel

# halt_on_error keeps the first report from being buried under later ones, and
# a sanitizer diagnostic must fail the run rather than be printed and ignored.
export ASAN_OPTIONS="halt_on_error=1:abort_on_error=1:detect_leaks=1"
export UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1"

ctest --test-dir "${build}" --output-on-failure
