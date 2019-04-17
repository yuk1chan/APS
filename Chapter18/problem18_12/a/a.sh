
g++ qwalk.cpp
./a.out < param.txt

sort psi.data -g > psi
rm psi.data
mv psi psi.data

gnuplot plot.plt

open psi.png

