#!/bin/sh
while read p;
do
  python3 TASK1.py --input "/geos/netdata/avtrain/data/3d/oosa/assignment/lvis/2009/$p" --output "$p" --res 100
done <2009.txt
