#!/bin/bash
set -xe
g++ -g lab4.cpp -o run -lm -lraylib
./run $1