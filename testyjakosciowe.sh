#!/bin/bash


input_dir="./testfiles"

input_files=("$input_dir"/*.txt)

for input_file in "${input_files[@]}"; do

    ./program.exe "$input_file"

done