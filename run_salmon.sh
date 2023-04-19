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

   # Run Salmon quant for the sample
    mkdir -p "salmon_output/$SAMPLE_NAME"
    $SALMON quant -i transcriptome_index -l A -1 $READ1 -2 $READ2 -o "salmon_output/$SAMPLE_NAME" --gcBias --validateMappings
done

# Aggregate the count data from each sample
find salmon_output -name "quant.sf" -exec tail -n +2 {} \; | sort -k1,1 | paste - - -d , > "temp_counts.txt"

# Calculate the number of samples
num_samples=$(find . -maxdepth 1 -type f -name "*trim.1.cor.fq.gz" | wc -l)

# Create the final count matrix file with the header
{
  echo "transcript_id$(printf ',%s' *trim.1.cor.fq.gz | sed 's/.trim.1.cor.fq.gz//g')"
  cat "temp_counts.txt" | awk -F, -v num_samples=$num_samples '{
    line = $1
    for (i = 1; i <= num_samples; i++) {
      line = line "," $(4 * i)
    }
    print line
  }'
} > $COUNT_MATRIX

# Remove the temporary file
rm "temp_counts.txt"







