#!/bin/sh

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -S . -B build

cmake --build build -j$(nproc)

sudo cmake --install build
