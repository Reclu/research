#!/bin/sh
echo "Building gif file..."
convert -delay 2 test.*.png anim.gif
echo "gif file built !"
echo "Building mp4 file..."
ffmpeg -f gif -i anim.gif video.mp4
echo "mp4 file built !"
