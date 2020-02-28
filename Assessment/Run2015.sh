#!/bin/sh
while read p;
do
  python3 TASK1.py --input "/geos/netdata/avtrain/data/3d/oosa/assignment/lvis/2015/$p" --output "$p.tif" --res 10
done <2015.txt
