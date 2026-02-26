#!/bin/env zsh

# Usage: ./compile.sh [build directory (default: ./build)]

if ! command -v cmake &> /dev/null; then
  echo "Missing CMake"
  exit 1
fi

if [[ $# < 1 ]]; then
  BUILD_DIR="build"
else
  BUILD_DIR=${1}
fi
SRC_DIR=$(pwd)/src


mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}
cmake ${SRC_DIR}
make
