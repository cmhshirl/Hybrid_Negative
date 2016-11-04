#!/bin/sh

for i in $(seq 1 30)
  do 
    time ./cell_lsodar $i &
    sleep 1
  done
  wait 
