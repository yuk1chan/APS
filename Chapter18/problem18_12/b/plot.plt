
psi(x) = exp(-x*x/2)/(sqrt(sqrt(pi)))
plot psi(x)

set terminal png
set output "psi.png"
replot "psi.data"

