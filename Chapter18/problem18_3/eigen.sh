
if [ -f "data.txt" ]; then
    rm data.txt
fi

if [ -f "phi_data.txt" ]; then
    rm phi_data.txt
fi

gcc eigen.c -o eigen

./eigen data.txt < input_data.txt

sort -n data.txt > phi_data.txt

gnuplot plot.plt
open phi.png


wc -l phi_data.txt > data_num.txt


if [ -f "data.txt" ]; then
    rm data.txt
fi


