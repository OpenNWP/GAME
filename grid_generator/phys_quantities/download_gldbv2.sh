#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# downloading GLDB data

wget "http://www.flake.igb-berlin.de/data/gldbv2.tar.gz"
tar -xzf gldbv2.tar.gz GlobalLakeDepth.dat
rm gldbv2.tar.gz
