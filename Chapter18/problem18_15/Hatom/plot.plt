set grid
set xrange[:10]
set yrange[0:0.5]
psi(x) = exp(-x)*x*x/sqrt(pi)
plot psi(x) title "psi"
replot "psi.data"
set terminal png
set output "psi.png"
replot
