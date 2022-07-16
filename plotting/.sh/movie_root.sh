#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

echo Creating movie ...
ffmpeg -y -hide_banner -loglevel warning -framerate 50 -i $fig_save_path/$run_id"_"$disp_shortname"_"$disp_level-%00ds.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p $fig_save_path/$run_id"_"$disp_shortname"_"$disp_level.mp4
if [ $? -ne 0 ]
then
  echo -e ${RED}Creaing movie failed.$NC
else
  echo Movie successfully created.
fi
