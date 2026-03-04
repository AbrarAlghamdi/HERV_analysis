#!/usr/bin/env bash

INDEX="/DataDrive/4TB/bladder_data/Cisplatin/salmon_index"

THREADS=16

for sample in */ ; do

    echo "Processing sample: $sample"

    cd "$sample"

    mkdir -p Results_salmon

    R1=$(ls *_R1_001.fastq.gz | sort | tr '\n' ' ')
    R2=$(ls *_R2_001.fastq.gz | sort | tr '\n' ' ')

    salmon quant \
        -i $INDEX \
        -l A \
        -1 $R1 \
        -2 $R2 \
        -p $THREADS \
        -o Results_salmon

    cd ..

done
