#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

if [ ! -d output ]
then
  mkdir output
fi

if [ ! -d nwp_init ]
then
  mkdir nwp_init
fi
