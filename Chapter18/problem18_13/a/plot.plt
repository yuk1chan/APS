set grid
set xrange[-5:5]
set yrange[0:2]
set terminal png
set output "psi.png"
plot "psi.data"
