#!/bin/sh

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -S . -B build

cmake --build build -j 6

sudo cmake --install build
