#!/bin/bash
set -xe
g++ -g lab3.cpp -o run -lm -lraylib
./run $1
