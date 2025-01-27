import pandas as pd


# Load SAM file into a DataFrame
def load_sam(filepath):
    """Loads a SAM file into a pandas DataFrame, only parsing alignment lines."""
    columns = [
        "QNAME",
        "FLAG",
        "RNAME",
        "POS",
        "MAPQ",
        "CIGAR",
        "RNEXT",
        "PNEXT",
        "TLEN",
        "SEQ",
        "QUAL",
        "TAGS",
    ]
    data = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith("@"):  # Skip headers
                fields = line.strip().split("\t")
                data.append(fields[:11] + [fields[11:]])  # Tags as a list
    return pd.DataFrame(data, columns=columns)


# Identify unique vs. chimeric alignments
def classify_alignments(df):
    """Classifies alignments as unique or chimeric based on TAGS and RNAME."""
    # Define TE categories based on RNAME
    unique_tes = [
        "Zebrafish_EnSpm-N49_DR#DNA/CACTA",
        "Zebrafish_EnSpm-N49B_DR#DNA/CACTA",
    ]  # Replace with actual names for unique TEs
    chimeric_te = (
        "Zebrafish-EnSpm-N49/N49B_DR"  # Replace with actual name for chimeric TE
    )

    # Define unique vs. chimeric based on RNAME
    df["Alignment_Type"] = df["RNAME"].apply(
        lambda x: "Unique"
        if x in unique_tes
        else "Chimeric"
        if x == chimeric_te
        else "Other"
    )

    # Check XA tag to detect if the read has multiple alignments, indicating potential chimeras
    df["Is_Chimera"] = df["TAGS"].apply(
        lambda tags: any("XA:Z:" in tag for tag in tags)
    )

    return df


# Count alignments per TE and chimeric proportions
def count_alignments(df):
    """Counts unique and chimeric alignments per TE and calculates chimeric proportions."""
    # Count total alignments by TE
    te_counts = df["RNAME"].value_counts().rename("Total_Alignments")

    # Count unique alignments by TE
    unique_counts = (
        df[df["Alignment_Type"] == "Unique"]["RNAME"]
        .value_counts()
        .rename("Unique_Alignments")
    )

    # Count chimeric alignments by TE
    chimeric_counts = (
        df[df["Is_Chimera"]]["RNAME"].value_counts().rename("Chimeric_Alignments")
    )

    # Combine into a summary DataFrame
    summary_df = pd.concat([te_counts, unique_counts, chimeric_counts], axis=1).fillna(
        0
    )
    summary_df["Chimeric_Proportion"] = (
        summary_df["Chimeric_Alignments"] / summary_df["Total_Alignments"]
    )

    return summary_df


# Analyze mismatches and gaps within CIGAR strings for further insights
def analyze_cigar(df):
    """Parses the CIGAR string to count mismatches and gaps for each alignment."""
    import re

    # Regular expressions to identify match (M), insertions (I), and deletions (D)
    match_re = re.compile(r"(\d+)M")
    insertion_re = re.compile(r"(\d+)I")
    deletion_re = re.compile(r"(\d+)D")

    # Functions to extract counts from the CIGAR string
    def count_matches(cigar):
        return sum(int(x) for x in match_re.findall(cigar))

    def count_insertions(cigar):
        return sum(int(x) for x in insertion_re.findall(cigar))

    def count_deletions(cigar):
        return sum(int(x) for x in deletion_re.findall(cigar))

    # Apply functions to extract values
    df["Matches"] = df["CIGAR"].apply(count_matches)
    df["Insertions"] = df["CIGAR"].apply(count_insertions)
    df["Deletions"] = df["CIGAR"].apply(count_deletions)

    return df


# Main analysis pipeline
def analyze_sam(filepath):
    """Executes the entire analysis pipeline on a SAM file."""
    # Step 1: Load and parse the SAM file
    df = load_sam(filepath)

    # Step 2: Classify alignments into unique or chimeric categories
    df = classify_alignments(df)

    # Step 3: Count unique and chimeric alignments per TE and calculate proportions
    summary_df = count_alignments(df)

    # Step 4: Analyze mismatches and gaps in CIGAR strings
    df = analyze_cigar(df)

    return df, summary_df


# Run the analysis and print results
sam_filepath = "../TE-piRNA_complementarity/piRNA_to_TE_filtered.sam"
df, summary_df = analyze_sam(sam_filepath)

# Display results
print("Alignment Summary by TE:")
print(summary_df)
print("\nDetailed Alignment Data (first few rows):")
print(df.head())


# Filter alignments for high-quality and high-identity hits
def filter_high_quality_alignments(df, mapq_threshold=30, min_match_ratio=0.9):
    """
    Filters alignments based on quality and match ratio.
    Args:
        df: DataFrame containing alignment data.
        mapq_threshold: Minimum MAPQ score to retain an alignment.
        min_match_ratio: Minimum ratio of matches to read length for high identity.
    """
    # Compute match ratio (matches / read length)
    df["Match_Ratio"] = df["Matches"] / df["SEQ"].str.len()

    # Apply filters
    filtered_df = df[
        (df["MAPQ"].astype(int) >= mapq_threshold)  # High MAPQ
        & (df["Match_Ratio"] >= min_match_ratio)  # High match identity
    ]

    return filtered_df


# Rewrite the table for clarity
def create_output_table(df):
    """
    Cleans and rewrites the alignment table for clarity and usability.
    Args:
        df: DataFrame containing alignment data.
    """
    # Map FLAG values to orientation
    df["Orientation"] = df["FLAG"].apply(
        lambda x: "Forward" if int(x) & 16 == 0 else "Reverse"
    )

    # Select and rename columns
    output_df = df[
        [
            "QNAME",
            "RNAME",
            "POS",
            "MAPQ",
            "Matches",
            "Insertions",
            "Deletions",
            "Match_Ratio",
            "Orientation",
        ]
    ].copy()

    output_df.rename(
        columns={
            "QNAME": "piRNA_Name",
            "RNAME": "Transposable_Element",
            "POS": "Start_Position",
            "MAPQ": "Mapping_Quality",
            "Matches": "Match_Count",
            "Insertions": "Insertion_Count",
            "Deletions": "Deletion_Count",
            "Match_Ratio": "Alignment_Identity",
        },
        inplace=True,
    )

    # Sort by Transposable Element and Alignment Identity
    output_df = output_df.sort_values(
        by=["Transposable_Element", "Alignment_Identity"], ascending=[True, False]
    )

    return output_df


# Apply the filtering and table creation
high_quality_df = filter_high_quality_alignments(df)
output_table = create_output_table(high_quality_df)

# Save the cleaned table to a file
output_table.to_csv(
    "../TE-piRNA_complementarity/High_Quality_TE-piRNA_Hits.csv", index=False
)

# Display the first few rows of the improved output
print(output_table.head())

import matplotlib.pyplot as plt
import seaborn as sns


def plot_unique_alignments_by_te(df):
    """Creates a bar plot showing unique alignment counts by TE, colored by sense/antisense orientation."""
    # Filter for unique alignments without mismatches across full piRNA
    unique_alignments = df[
        (df["Alignment_Type"] == "Unique")
        & (df["Insertions"] == 0)
        & (df["Deletions"] == 0)
    ]
    unique_alignments["Orientation"] = unique_alignments["FLAG"].apply(
        lambda x: "Sense" if int(x) & 16 == 0 else "Antisense"
    )

    # Count unique alignments per TE and orientation
    alignment_counts = (
        unique_alignments.groupby(["RNAME", "Orientation"])
        .size()
        .reset_index(name="Counts")
    )

    # Plot
    plt.figure(figsize=(10, 6))
    plot = sns.barplot(
        data=alignment_counts,
        x="RNAME",
        y="Counts",
        hue="Orientation",
        palette="viridis",
    )
    plot.set_yscale("log")
    _ = plot.set(xlabel="Class", ylabel="Survived")
    plt.xticks(rotation=45, ha="right")
    plt.title("Unique Alignments by TE and Orientation (Sense/Antisense)")
    plt.xlabel("Target TE")
    plt.ylabel("log(Unique Alignment Counts)")
    plt.legend(title="Orientation")
    plt.tight_layout()
    plt.show()


plot_unique_alignments_by_te(df)
