#!/bin/bash

input_dir="./errors"

input_files=("$input_dir"/*.txt)

for input_file in "${input_files[@]}"; do
    echo "Przetwarzanie pliku: $input_file"
    
    while true; do
        ./program.exe "$input_file"

        # Uruchomienie drugiego programu i sprawdzenie komunikatu OK
        ./ckrptw.exe "$input_file" wyniki.txt


    done
done
