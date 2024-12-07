#!/bin/bash

input_dir="./testfiles"

input_files=("$input_dir"/*.txt)

OK=0
NO_OK=0
ALL=0
NO_OK_FILES=() # Tablica plików z wynikami NO_OK

for i in $(seq 1);do
    for input_file in "${input_files[@]}"; do
        echo "Przetwarzanie pliku: $input_file"
        ((ALL++))
        ./program.exe "$input_file"
        output_program=$(./ckrptw.exe "$input_file" wyniki.txt)

        if [[ "$output_program" == *OK* ]]; then
            ((OK++))
        else
            ((NO_OK++))
            NO_OK_FILES+=("$input_file") # Dodanie pliku do tablicy
        fi
    done
done

echo "OK: $OK"
echo "NO_OK: $NO_OK"

if ((NO_OK > 0)); then
    echo "Pliki z wynikami NO_OK:"
    for file in "${NO_OK_FILES[@]}"; do
        echo "$file"
    done
else
    echo "Wszystkie pliki przetworzone poprawnie."
fi
