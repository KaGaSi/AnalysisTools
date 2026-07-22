#!/bin/env zsh

cmake -S src -B build-san -DSANITIZE=ON

cmake --build build-san

ASAN_OPTIONS=abort_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
  ctest --test-dir build-san --output-on-failure
