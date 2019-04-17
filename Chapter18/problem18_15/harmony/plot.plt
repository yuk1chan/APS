set grid
set title "V(x) = x^2/2 + 0.00000x^3"
set xrange[-5:5]
set yrange[0:2]
V(x) = x*x/2 + 0.00000*x*x*x
plot V(x)
set terminal png
set output "psi.png"
replot "psi.data"
