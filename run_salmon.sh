#!/bin/bash

# Path to the transcriptome fasta file
TRANSCRIPTOME="path/to/transcriptome.fasta"

# Path to the Salmon executable
SALMON="/path/to/salmon"

# Index the transcriptome
$SALMON index -t $TRANSCRIPTOME -i transcriptome_index

# Initialize the output files for the count matrix and sample metadata
COUNT_MATRIX="Count_Matrix.csv"
SAMPLE_METADATA="Sample_Metadata.csv"

# Write the header for the count matrix and sample metadata
echo "transcript_id$(printf ',%s' *trim.1.cor.fq.gz | sed 's/.trim.1.cor.fq.gz//g')" > $COUNT_MATRIX
echo "sample,condition" > $SAMPLE_METADATA

# Function to classify samples based on their filenames
classify_condition() {
    if [[ $1 == *"[Ss]aline"* && $1 == *"[Ff]ooddep"* ]]; then
        echo "saline_fooddep"
    elif [[ $1 == *"[Ss]aline"* && ($1 == *"[Aa]d"*"lib"* || $1 == *"[Aa]dlib"*) ]]; then
        echo "saline_adlib"
    elif [[ $1 == *"[Ll]eptin"* && $1 == *"[Ff]ooddep"* ]]; then
        echo "leptin_fooddep"
    elif [[ $1 == *"[Ll]eptin"* && ($1 == *"[Aa]d"*"lib"* || $1 == *"[Aa]dlib"*) ]]; then
        echo "leptin_adlib"
    else
        echo "Unknown"
    fi
}

# Iterate over the files and run salmon quant for each paired-read sample
for READ1 in *trim.1.cor.fq.gz; do
    # Find the corresponding file with .2
    READ2=${READ1/trim.1.cor.fq.gz/trim.2.cor.fq.gz}

    # Extract the sample name and condition
    SAMPLE_NAME=$(echo $READ1 | cut -d '.' -f1)
    CONDITION=$(classify_condition $SAMPLE_NAME)

    # Add sample metadata to the sample metadata file
    echo "${SAMPLE_NAME},${CONDITION}" >> $SAMPLE_METADATA

    # Create a directory for the output of the sample
    mkdir -p $SAMPLE_NAME

    # Run Salmon quant for the sample
    $SALMON quant -i transcriptome_index -l A -1 $READ1 -2 $READ2 -o $SAMPLE_NAME --gcBias --validateMappings
done

# Function to extract count from a single file
extract_count() {
    grep -m1 -w $1 $2/quant.sf | awk '{print $4}'
}

# Aggregate the count data from each sample
for TRANSCRIPT in $(cut -f1 ${SAMPLE_NAME}/quant.sf | tail -n +2); do
    COUNTS=$(paste -d, <(printf "${TRANSCRIPT}") <(paste -d, -s <(for SAMPLE_DIR in *trim.1.cor.fq.gz; do
        SAMPLE_NAME=$(echo $SAMPLE_DIR | cut -d '.' -f1)
        extract_count $TRANSCRIPT $SAMPLE_NAME
    done)))
    echo "${COUNTS}" >> $COUNT_MATRIX
done
