#!/bin/bash

# File paths
PIRNA_FILE="drer_pirna.fa"
TE_FILE="transposon_seq.fa"
OUTPUT_SAM="piRNA_to_TE_alignment.sam"
OUTPUT_BAM="piRNA_to_TE_alignment.bam"
FILTERED_BAM="piRNA_to_TE_filtered.bam"

# Step 1: Index the TE sequences
echo "Indexing transposable element sequences..."
bwa index "$TE_FILE"

# Step 2: Align piRNAs to TE sequences
echo "Aligning piRNA sequences to transposable element sequences..."
bwa aln "$TE_FILE" "$PIRNA_FILE" > piRNA_to_TE.sai
bwa samse "$TE_FILE" piRNA_to_TE.sai "$PIRNA_FILE" > "$OUTPUT_SAM"

# Step 3: Convert SAM to BAM, sort, and filter for high-quality alignments
echo "Converting SAM to sorted BAM..."
samtools view -Sb "$OUTPUT_SAM" | samtools sort -o "$OUTPUT_BAM"

echo "Filtering for high-quality alignments (MAPQ >= 20)..."
samtools view -h -b -q 20 "$OUTPUT_BAM" > "$FILTERED_BAM"

# Step 4: Index the filtered BAM file for viewing
echo "Indexing the filtered BAM file..."
samtools index "$FILTERED_BAM"

# Clean up intermediate files
echo "Cleaning up intermediate files..."
rm piRNA_to_TE.sai "$OUTPUT_SAM"

echo "Alignment and filtering complete. Results in $FILTERED_BAM and $FILTERED_BAM.bai"
