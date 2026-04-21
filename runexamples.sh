#!/bin/bash
# Run all examples
./build/bin/ex1
./build/bin/ex2
./build/bin/ex3
./build/bin/ex7

# go to the plotting
cd plotting/functions
python3 ex1.py &
python3 ex2.py &
python3 ex3.py &
python3 ex7.py &
cd ../../