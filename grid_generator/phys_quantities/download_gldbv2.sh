#!/bin/bash

# downloading GLDB data

cd phys_quantities
wget "http://www.flake.igb-berlin.de/data/gldbv2.tar.gz"
tar -xzf gldbv2.tar.gz GlobalLakeStatus.dat
rm gldbv2.tar.gz
cd ..
