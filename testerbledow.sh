#!/bin/bash

input_dir="./testfiles"

input_files=("$input_dir"/*.txt)

while true; do
    for input_file in "${input_files[@]}"; do
        echo "Przetwarzanie pliku: $input_file"
    
        ./program.exe "$input_file"

        # Uruchomienie drugiego programu i sprawdzenie komunikatu OK
        ./ckrptw.exe "$input_file" wyniki.txt


    done
done
