
gcc normalization.c -o normalization

if [ -f "normalization.txt" ]; then
    rm normalization.txt
fi

./normalization data_num.txt phi_data.txt normalization.txt


