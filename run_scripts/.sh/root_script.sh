#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

run_dir=$game_home_dir/output/$run_id
if [ -d $run_dir ]
then
  rm -r $run_dir
fi
mkdir $run_dir
cd $run_dir

if [ ! -f $game_home_dir/build/game ]
then
  echo "Executable game missing. Compile first. Aborting run."
  cd - > /dev/null
  exit 1
fi

$game_home_dir/build/game

cd - > /dev/null
