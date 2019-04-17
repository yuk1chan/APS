set grid
set title "V(x) = x^2/2 + 0.10000x^3"
set xrange[-10:10]
set yrange[0:2]
V(x) = x*x/2 + 0.10000*x*x*x
plot V(x)
set terminal png
set output "psi.png"
replot "psi.data"
