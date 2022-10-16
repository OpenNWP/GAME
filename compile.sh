#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

if [ ! -d build ]
then
  mkdir build
fi

cd build

d_value=False
f_value=False
c_value=False
while getopts "dfc" opt; do
  case $opt in
    d) # debugging flag
      d_value=True
      ;;
    f) # aggressive optimization flag
      f_value=True
      ;;
    c) # compile-time configuration flag
      c_value=True
      ;;
    \?)
      echo "Invalid option: -$OPTARG. Compiling anyway."
      ;;
  esac
done

cmake -DDEBUGGING=$d_value -DFAST=$f_value -DCOMPILE_TIME_CONFIG=$c_value ..
make

cd ..
