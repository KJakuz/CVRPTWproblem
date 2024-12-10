#!/bin/bash


input_dir="./testyjakosci"

input_files=("$input_dir"/*.txt)

for i in $(seq 1 2 3 4 5);do
    for input_file in "${input_files[@]}"; do

        ./program.exe "$input_file"

    done
done