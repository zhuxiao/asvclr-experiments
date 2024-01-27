#!/bin/bash

input_file="SRR_Acc_List.txt"
output_file="SRR885_whole.fastq"
> "$output_file"

if [ ! -f "$input_file" ]; then
    echo "Input file not found."
    exit 1
fi
mkdir SRX5327410
while IFS= read -r line;do
    cd SRX5327410
    echo "$line downloading..." 
    prefetch $line
    echo "$line download successfully"
    cd $line
    echo "fasterq-dump $line"
    fasterq-dump $line.sra
    echo "done."
    echo "cat starting"
    cat $line.fastq >> ../../$output_file
    echo "cat done"
    cd ../../
done < "$input_file"
echo "All file are download and dump successfully"


