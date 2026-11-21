#!/usr/bin/env bash

set -euo pipefail

BASE_URL="https://zenodo.org/records/800651/files"

mkdir -p data/raw data/reference

samples=(
    F3D0 F3D1 F3D2 F3D3 F3D5
    F3D141 F3D142 F3D143 F3D144 F3D145
    Mock
)

for sample in "${samples[@]}"; do
    for read in R1 R2; do

        filename="${sample}_${read}.fastq"

        if [[ ! -s "data/raw/$filename" ]]; then
            echo "Downloading $filename"

            curl -fL --retry 3 \
                "$BASE_URL/$filename?download=1" \
                -o "data/raw/$filename"
        else
            echo "Already downloaded: $filename"
        fi

    done
done

references=(
    silva.v4.fasta
    trainset9_032012.pds.fasta
    trainset9_032012.pds.tax
)

for filename in "${references[@]}"; do

    if [[ ! -s "data/reference/$filename" ]]; then
        echo "Downloading $filename"

        curl -fL --retry 3 \
            "$BASE_URL/$filename?download=1" \
            -o "data/reference/$filename"
    else
        echo "Already downloaded: $filename"
    fi

done

echo "All downloads completed."
