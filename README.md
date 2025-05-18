# This repository contains all analyses performed on *D. rerio* *ich* mutant strain

Screening for transposable element insertions genome-wide was performed using BLAST search algorithm:

```
# Prepare genome database
makeblastdb -in [strain_genome.fa] -dbtype nucl -out [strain_genome_db]
```

```
# Perform search
blastn  -query fasta/transposon_seq.fa \
        -db [strain_genome_db] \
        -out [strain name]_crossreff.tsv \
        -outfmt 6
```

The structure of the repository is the following:

- Screening for transposable element (TE) insertion using BLAST (outputs in `strain_hits/`)
  - Multiple strains were screened for the N49, N49B, and N49/N49B insertion
  - For plotting purposes the chromosome names from each strain had to be crossreferenced (`[strain name]_crossreff.tsv`)
  - Output of TE insertion counts per 2 Mb bins are also stored here (`[strain name]-TE_hits_binned.csv`)  
- Different plots can be found in the `outputs/`
- Fasta files used during analyses are stored in `fasta/`
- Scripts used:
  - `scripts/fish_transposon_analysis.R` - mains script used for plotting and analysis
  - `scripts/2025Varga-ich_mapping_ms-figure_scripts.R` - script used for generating publication plots


