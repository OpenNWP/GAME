#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

d_value=False
f_value=False
c_value=False
s_value=False
while getopts "dfcs" opt; do
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
    s) # single-precision flag
      s_value=True
      ;;
    \?)
      echo "Invalid option: -$OPTARG. Compiling anyway."
      ;;
  esac
done

if [ ! -d build ]
then
  mkdir build
fi

cd build

cmake -DDEBUGGING=$d_value -DFAST=$f_value -DCOMPILE_TIME_CONFIG=$c_value -DSINGLE_PRECISION=$s_value ..
make

cd ..
