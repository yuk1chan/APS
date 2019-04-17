set terminal png
set grid

set out "phi.png"
plot "phi_data.txt" using 1:2 w lp

