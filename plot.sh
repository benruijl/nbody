#!/bin/bash

gcc nbody.c -lm -o nbody
./nbody
gnuplot plot.conf
