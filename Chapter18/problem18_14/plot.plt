set grid
set title "V(x) = 5.000000"
set xrange[-10:10]
set yrange[0:2]
set terminal png
set output "psi.png"
plot "psi.data"
