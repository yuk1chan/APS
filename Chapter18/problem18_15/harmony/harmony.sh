#!/bin/bash

g++ GFMC.cpp
./a.out < param.txt

gnuplot plot.plt

open psi.png


