#!/bin/bash

gcc nbody.c -lm -o nbody
./nbody nbody.conf
gnuplot -persist plot.conf
