#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

d_value=False
while getopts "d" opt; do
  case $opt in
    d) # debugging flag
      d_value=True
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

cmake -DDEBUGGING=$d_value ..
make

cd ..
