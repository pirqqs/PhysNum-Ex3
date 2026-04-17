#!/bin/bash

set -e


echo "Compiling C++ code..."
g++ -std=c++17 -O2 main.cc -o main

echo "Running parameter scan..."
python3 ParameterScan.py

echo "Generating figures..."
python3 figures.py

echo "Done."