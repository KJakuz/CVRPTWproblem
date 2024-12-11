#!/bin/bash


input_dir="./testyjakosci"

input_files=("$input_dir"/*.txt)

for i in $(seq 1 10);do
    for input_file in "${input_files[@]}"; do

        ./program150.exe "$input_file"

    done
done