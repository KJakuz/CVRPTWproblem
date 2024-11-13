#!/bin/bash


input_dir="./testfiles"

input_files=("$input_dir"/*.txt)

for input_file in "${input_files[@]}"; do
    echo "Przetwarzanie pliku: $input_file"

    ./program.exe "$input_file"
    ./ckrptw.exe "$input_file" wyniki.txt

done
