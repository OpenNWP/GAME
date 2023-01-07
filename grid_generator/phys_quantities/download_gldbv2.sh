#!/bin/bash

# downloading GLDB data

wget "http://www.flake.igb-berlin.de/data/gldbv2.tar.gz"
tar -xzf gldbv2.tar.gz GlobalLakeDepth.dat
rm gldbv2.tar.gz
