
if [ -f "data.txt" ]; then
    rm data.txt
fi

if [ -f "phi0_data.txt" ]; then
    rm phi0_data.txt
fi

# コンパイル
gcc calc_phi0.c -o phi0
gcc calc_phib.c -o phib
gcc calc_integral.c -o integral

# phi0
./phi0 data.txt < input_phi0_data.txt

sort -n data.txt > phi0_data.txt

# phibのためにdata.txtを削除
if [ -f "data.txt" ]; then
    rm data.txt
fi

# phib
./phib data.txt < input_phib_data.txt

sort -n data.txt > phib_data.txt


#gnuplot plot.plt
#open phi.png

wc -l phi0_data.txt > data_phi0_num.txt
wc -l phib_data.txt > data_phib_num.txt

./integral data_phi0_num.txt phi0_data.txt data_phib_num.txt phib_data.txt


if [ -f "data.txt" ]; then
    rm data.txt
fi

