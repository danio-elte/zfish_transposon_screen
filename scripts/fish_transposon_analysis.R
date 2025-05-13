# Load required libraries for analysis and visualization
library(GenomicRanges)                  
library(rtracklayer)                        
library(RColorBrewer)                       
library(tidyverse)                         
library(Biostrings)                         
library(patchwork)                          
library(circlize)                           
library(biomaRt)                            
library(ChIPseeker)                         
library(TxDb.Drerio.UCSC.danRer11.refGene)
library(flextable)
library(officer)
library(officedown) 
library(magrittr)
library(universalmotif)
library(seqLogo)
library(tidyplots)

# Set working directory for relative paths
setwd(this.path::here())

# Load the genome annotation database for zebrafish
# This database provides transcript and genomic annotations.
txdb <- TxDb.Drerio.UCSC.danRer11.refGene

# Read BLAST output file (-outfmt 6) and transposon sequences
# Filters results based on e-value threshold (≤ 1e-10)
blastout <- read_tsv("../strain_hits/Tuebingen.out", col_names = FALSE) %>%
  filter(X11 <= 1e-10)

# Load cross-reference file for chromosome mapping
chr_crossref <- read_tsv("../strain_hits/Tuebingen_crossreff.tsv")

# Read transposon sequences from a FASTA file to extract lengths
seqs <- readDNAStringSet("../fasta/transposon_seq.fa")

# Create a data frame containing transposon names and their lengths
expect <- data.frame(name = names(seqs), length = width(seqs))

# Initialize an empty data frame to store hit proportions for analysis
hits <- data.frame(percent = as.double(),
                   hybrid = as.double(),
                   N49B = as.double(),
                   N49 = as.double())

# Define proportion thresholds for transposon alignment
prop <- seq(0, 1, 0.05)

# Loop through each proportion threshold to calculate transposon hits
for (i in 1:length(prop)) {
  prop_temp <- prop[i]
  
  # Adjust the threshold for the first iteration
  if (i == 1) prop_temp <- prop_temp + 0.01

  # Calculate the proportion of aligned transposons
  dims <- blastout %>%
    left_join(expect, by = c("X1" = "name")) %>%
    mutate(proportion = (X8 - X7) / length) %>%
    filter(proportion >= prop_temp) %>%
    group_by(X1) %>%
    summarise(N = n())

  # Store the proportion and hits for specific transposon types
  hits[i, 1] <- prop_temp
  hits[i, 2] <- ifelse("Zebrafish-EnSpm-N49/N49B_DR" %in% dims$X1, 
                       dims %>% filter(X1 == "Zebrafish-EnSpm-N49/N49B_DR") %>% pull(N), 0)
  hits[i, 3] <- ifelse("Zebrafish_EnSpm-N49B_DR#DNA/CACTA" %in% dims$X1, 
                       dims %>% filter(X1 == "Zebrafish_EnSpm-N49B_DR#DNA/CACTA") %>% pull(N), 0)
  hits[i, 4] <- ifelse("Zebrafish_EnSpm-N49_DR#DNA/CACTA" %in% dims$X1, 
                       dims %>% filter(X1 == "Zebrafish_EnSpm-N49_DR#DNA/CACTA") %>% pull(N), 0)
}

# Replace any NA values in the hits data frame with zeros
hits[is.na(hits)] <- 0

# Plot the number of transposon hits across different proportion thresholds
p1 <- hits %>%
  tibble() %>%
  pivot_longer(cols = 2:4) %>%
  mutate(name = case_when(name == "N49" ~ "EnSpm-N49_DR",
                          name == "N49B" ~ "EnSpm-N49B_DR",
                          name == "hybrid" ~ "EnSpm-N49/N49B_DR"),
         name = factor(name, levels = c("EnSpm-N49_DR", "EnSpm-N49B_DR", "EnSpm-N49/N49B_DR")),
         percent = percent * 100) %>%
  ggplot(aes(x = percent, y = value, color = name)) +
  geom_line(linewidth = 4) +
  xlab("% of transposon aligned to genome") +
  ylab("Number of transposon hits") +
  scale_color_manual(values = c("#FFAC94", "#FF83B1", "#7FAF75")) +
  labs(color = NULL) +
  geom_vline(xintercept = 80, linetype = "dashed", linewidth = 2) +
  theme_linedraw() +
  theme(
    legend.position = c(0.78, 0.86),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    text = element_text(size = 20),
    axis.text = element_text(size = 18)
  )


blastout %>%
  left_join(expect, by = c("X1" = "name")) %>%
  mutate(proportion = (X8 - X7) / length) %>%
  filter(proportion >= 0.8) %>%
  mutate(X2 = paste0("chr", chr_crossref$`Chromosome name`[match(.$X2, chr_crossref$`RefSeq seq accession`)])) %>%
  group_by(X1) %>%
  summarise(N = n())

# Plot BLAST hits at 80% alignment threshold grouped by chromosome
p2 <- blastout %>%
  left_join(expect, by = c("X1" = "name")) %>%
  mutate(proportion = (X8 - X7) / length) %>%
  filter(proportion >= 0.8) %>%
  mutate(X2 = paste0("chr", chr_crossref$`Chromosome name`[match(.$X2, chr_crossref$`RefSeq seq accession`)])) %>%
  group_by(X1, X2) %>%
  summarise(N = n()) %>%
  mutate(X2 = as.double(str_remove_all(X2, "chr"))) %>%
  arrange(X2) %>%
  mutate(X1 = case_when(X1 == "Zebrafish_EnSpm-N49_DR#DNA/CACTA" ~ "EnSpm-N49_DR",
                          X1 == "Zebrafish_EnSpm-N49B_DR#DNA/CACTA" ~ "EnSpm-N49B_DR",
                          X1 == "Zebrafish-EnSpm-N49/N49B_DR" ~ "EnSpm-N49/N49B_DR"),
         X1 = factor(X1, levels = c("EnSpm-N49_DR", "EnSpm-N49B_DR", "EnSpm-N49/N49B_DR"))) %>%
  ggplot(aes(x = as.factor(X2), y = N, fill = X1)) +
  geom_bar(stat = 'identity', color = 'black') +
  facet_wrap(~X1, ncol = 1) +
  theme_linedraw() +
  xlab("Chromosome") +
  ylab("Blast hits with 80% of\ntransposon aligned to genome") +
  scale_fill_manual(values = c("#FFAC94", "#FF83B1", "#7FAF75")) +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(colour = 'black'),
    text = element_text(size = 15)
  )

# Combine and save the two plots
p3 <- p1 + p2 + plot_annotation(tag_levels = "A")
ggsave("../output/blast_summary.png", p3, width = 14, height = 10, units = 'in', dpi = 300)

# Save p1
ggsave("../output/transposon_hits.pdf", p1, width = 7, height = 5, units = 'in', dpi = 300)

##########################
### Strain Comparisons ###
##########################

# Define file paths
hit_files <- list.files("../strain_hits/", pattern = "\\.out$", full.names = TRUE)
crossref_files <- list.files("../strain_hits/", pattern = "crossreff", full.names = TRUE)

# Initialize hits data frame
hits <- data.frame()

# Loop through files and process hits
for (i in seq_along(hit_files)) {
    # Read and process the hit files
    temp_hits <- read_tsv(hit_files[i], col_names = FALSE, skip_empty_rows = TRUE)
    colnames(temp_hits) <- c(
        "query acc.ver", "subject acc.ver", "% identity", "alignment length",
        "mismatches", "gap opens", "q. start", "q. end",
        "s. start", "s. end", "evalue", "bit score"
    )
    
    # Ensure all necessary columns are character
    temp_hits <- temp_hits %>%
        mutate(across("query acc.ver":"subject acc.ver", as.character))
    
    # Process cross-reference files
    temp_crossref <- read_tsv(crossref_files[i], col_names = TRUE, skip_empty_rows = TRUE) %>%
        mutate(across(everything(), as.character))
    
    # Perform join and append
    if (str_detect(hit_files[i], "Tuebingen")) {
        temp_hits <- full_join(temp_hits, temp_crossref, by = c("subject acc.ver" = "RefSeq seq accession"))
        strain_name <- str_remove_all(basename(hit_files[i]), "\\.out$")
        temp_hits <- temp_hits %>%
            mutate(strain = strain_name) %>%
            dplyr::select(-`GenBank seq accession`)
    } else {
        temp_hits <- full_join(temp_hits, temp_crossref, by = c("subject acc.ver" = "GenBank seq accession"))
        strain_name <- str_remove_all(basename(hit_files[i]), "\\.out$")
        temp_hits <- temp_hits %>%
            mutate(strain = strain_name) %>%
            dplyr::select(-`RefSeq seq accession`)
    }
    
    # Append to the main hits data frame
    hits <- bind_rows(hits, temp_hits)
}


# Clean and process hits
hits <- hits %>%
  filter(!is.na(`s. start`), !is.na(`s. end`)) %>%
  mutate(
    `s. start` = pmin(`s. start`, `s. end`),
    `s. end` = pmax(`s. start`, `s. end`)
  ) %>%
  left_join(expect, by = c("query acc.ver" = "name")) %>%
  filter(
    evalue <= 1e-10,
    (`q. end` - `q. start`) / `length` >= 0.8,
    `Chromosome name` != "NA"
  ) %>%
  mutate(
    `Chromosome name` = factor(`Chromosome name`, levels = c(seq(1, 25), "Un")),
    `query acc.ver` = recode(
      `query acc.ver`,
      "Zebrafish-EnSpm-N49/N49B_DR" = "EnSpm-N49/N49B",
      "Zebrafish_EnSpm-N49B_DR#DNA/CACTA" = "EnSpm-N49B",
      "Zebrafish_EnSpm-N49_DR#DNA/CACTA" = "EnSpm-N49"
    )
  )

# Plot hits by chromosome
p <- hits %>%
  mutate(`query acc.ver` = factor(`query acc.ver`, levels = c("EnSpm-N49", "EnSpm-N49B", "EnSpm-N49/N49B"))) %>%
  ggplot(aes(x = `Chromosome name`, fill = `query acc.ver`)) +
  geom_bar(stat = "count", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("#FFAC94", "#FF83B1", "#7FAF75")) +
  facet_grid(strain ~ `query acc.ver`) +
  theme_linedraw() +
  xlab("Chromosome") +
  ylab("Quality-filtered BLAST\nhit counts") +
  theme(
    legend.position = "none",
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(colour = "black"),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  )
ggsave("output/strain_hits.png", p, width = 9, height = 5, units = "in", dpi = 300)

# Prepare GRanges object
hits_gr <- GRanges(
  seqnames = paste0("chr", hits$`Chromosome name`),
  ranges = IRanges(start = hits$`s. start`, end = hits$`s. end`),
  strain = hits$strain
)

ctnnb2 <- GRanges(seqnames = "chr19", IRanges(start = 47975912, end = 47997421))

# Find overlaps with CTNNB2 gene  
hits_ctnnb2 <- findOverlaps(hits_gr, ctnnb2)
hits_ctnnb2 <- hits[queryHits(hits_ctnnb2), ]
hits_ctnnb2

####################
### Circos Plots ###
####################

# Prefix chromosome names
hits$`Chromosome name` <- str_replace_all(hits$`Chromosome name`, "^", "chr")
hits <- hits %>% 
  filter(str_detect(`Chromosome name`, "^chr[0-9]+$"))

# Refactor Chromosome name levels to start from chr1 and end at chr25
hits$`Chromosome name` <- factor(
  hits$`Chromosome name`,
  levels = paste0("chr", 1:25),
  ordered = TRUE
)

# Filter to include only chromosomes within the range chr1 to chr25
hits <- hits %>%
  filter(!is.na(`Chromosome name`))

# Define bin size
bin_size <- 2e6  # 2Mb bins

# Process hits for each strain separately
strain_list <- split(hits, hits$strain)
strain_results <- list()

for (strain in names(strain_list)) {
  strain_hits <- strain_list[[strain]]
  
  genome_inf <- read_tsv(paste0("strain_hits/", strain, "_crossreff.tsv")) %>%
    mutate(`Chromosome name` = paste0("chr", `Chromosome name`),
           Start = 1,
           End = `Seq length`) %>%
    filter(str_detect(`Chromosome name`, "chr[0-9]+")) %>%
    dplyr::select(`Chromosome name`, Start, End)
    
  
  # Create chromosome-specific bins based on strain-specific chromosome lengths
  genome_bins <- genome_inf %>%
    group_by(`Chromosome name`) %>%
    summarize(
      Start = seq(1, max(as.integer(End)), by = bin_size),
      End = c(seq(bin_size, max(as.integer(End)), by = bin_size), max(as.integer(End)))
    ) %>%
    ungroup() %>%
    mutate(Bin = row_number())
  
  # Convert strain_hits and genome_bins to GRanges
  hit_gr <- GRanges(
    seqnames = strain_hits$`Chromosome name`,
    ranges = IRanges(start = strain_hits$`s. start`, end = strain_hits$`s. end`),
    TE = strain_hits$`query acc.ver`
  )
  
  bin_gr <- GRanges(
    seqnames = genome_bins$`Chromosome name`,
    ranges = IRanges(start = genome_bins$Start, end = genome_bins$End),
    Bin = genome_bins$Bin
  )
  
  # Map hits to bins
  hits_to_bins <- findOverlaps(hit_gr, bin_gr)
  hits_mapped <- data.frame(
    Bin = mcols(bin_gr)$Bin[subjectHits(hits_to_bins)],
    TE = mcols(hit_gr)$TE[queryHits(hits_to_bins)]
  )
  
  # Count hits by TE and Bin
  te_bins <- hits_mapped %>%
    group_by(TE, Bin) %>%
    summarize(hit_count = n(), .groups = "drop") %>%
    right_join(genome_bins, by = "Bin") %>%
    mutate(hit_count = replace_na(hit_count, 0))
  
  # Store the result for this strain
  strain_results[[strain]] <- list(
    genome_bins = genome_bins,
    te_bins = te_bins
  )
}

# Generate circos plots for each strain
for (strain in names(strain_results)) {
  strain_data <- strain_results[[strain]]
  
  genome_bins <- strain_data$genome_bins
  genome_bins$`Chromosome name` <- factor(
    genome_bins$`Chromosome name`,
    levels = paste0("chr", 1:25),
    ordered = TRUE
  )
  te_bins <- strain_data$te_bins %>%
    filter(!(is.na(TE))) %>%
    mutate(TE = factor(TE, levels = c("EnSpm-N49", "EnSpm-N49B", "EnSpm-N49/N49B")))
  
  te_bins <- genome_bins %>%
    distinct(`Chromosome name`, Start, End, Bin) %>%  # Get unique chromosome bins
    expand(
      `Chromosome name`, Start, Bin, TE=unique(te_bins$TE)) %>% 
    mutate(End = Start + 2e6 - 1) %>% 
    left_join(te_bins, by = c("Chromosome name", "Start", "End", "Bin", "TE")) %>% 
    mutate(hit_count = ifelse(is.na(hit_count), 0, hit_count))
  
  # Prepare colors for TEs
  te_colors <- setNames(c("#FFAC94", "#FF83B1", "#7FAF75"), 
                        levels(te_bins$TE))
  
  # Define factors and xlim
  chromosome_levels <- unique(genome_bins$`Chromosome name`)
  xlim <- genome_bins %>%
    group_by(`Chromosome name`) %>%
    summarize(Start = min(Start), End = max(End)) %>%
    ungroup() %>%
    dplyr::select(Start, End) %>%
    as.matrix()
  
  # Export DF
  write.csv(te_bins %>% dplyr::select(-Bin) %>% distinct(.keep_all = T), quote = F, row.names = F, 
            file = paste0('strain_hits/', strain, '_TE-hits_binned.csv'))
  
  
  # Create the circos plot
  pdf(paste0("output/circos_", strain, ".pdf"), width = 9, height = 9)
  circos.clear()

  circos.par(
    start.degree = 90, 
    gap.degree = c(rep(1, length(chromosome_levels) - 1), 3),
    points.overflow.warning = FALSE)

  circos.initialize(
    factors = chromosome_levels,
    xlim = cbind(xlim[, 1], xlim[, 2])
  )

  
  # Add tracks for each TE
  for (te in levels(te_bins$TE)) {
    te_data <- te_bins %>% 
      filter(TE == !!te) %>%
      dplyr::select(-Bin) %>%
      distinct(.keep_all = T)
    
    circos.trackPlotRegion(
      factors = te_data$`Chromosome name`,
      y = te_data$hit_count,
      x = te_data$Start,
      track.height = 0.1,
      ylim = c(0, 4),
      panel.fun = function(region_start, y) {
                if(CELL_META$track.index == 1){
          circos.xaxis(
            #h = 'bottom',
            labels.facing = 'clockwise',
            labels.niceFacing = TRUE,
            major.at = c(0, 20e6, 40e6, 80e6),  # Positions in base pairs
            major.tick.length = 0.8,  # Adjust the tick size
            labels = c("0 Mb", "20 Mb", "40 Mb", "80 Mb"),  # Labels in Mb
            labels.cex = 0.6,  # Font size for labels
            sector.index = CELL_META$sector.index,  # Automatically match sector
          )
        }
        circos.barplot(y, region_start, col = te_colors[[te]], border = te_colors[[te]]) 
      }
        )

  }
  
  # Add labels
  for (chromosome in unique(genome_bins$`Chromosome name`)) {
    mid_point <- genome_bins %>%
      filter(`Chromosome name` == chromosome) %>%
      summarize(mid = mean(c(Start, End)))
    
    circos.text(
      sector.index = chromosome,
      track.index = 3,
      x = mid_point$mid,
      y = max(te_bins$hit_count) * -0.9,  # Adjust y for placement
      labels = chromosome,
      cex = 0.8,
      facing = "bending.inside"
    )
  }
  
  text(0, 0, paste0("Strain:\n", strain), cex = 1.2)
  legend(x = -0.23, y = -0.2,
         pch = 16, 
         legend = names(te_colors), 
         col = te_colors, 
         pt.cex = 1.5, 
         cex = 1, 
         text.col = "black", 
         horiz = F , 
         inset = c(0.1, 0.1))

# Y-axis labels
circos.yaxis(sector.index = "chr1",
               at = c(seq(1, 5, 1)), 
               track.index = 1,
               labels.cex = 0.5)
circos.yaxis(sector.index = "chr1",
               at = c(seq(1, 5, 1)), 
               track.index = 2,
               labels.cex = 0.5)
circos.yaxis(sector.index = "chr1",
               at = c(seq(1, 5, 1)), 
               track.index = 3,
               labels.cex = 0.5)




  circos.clear()
  dev.off()
}


###################################
### RefGenome with Gene Density ###
###################################

# Reuse the existing data loading code
ensembl <- useEnsembl(biomart = "genes",
                      dataset = "drerio_gene_ensembl")

genes <- getBM(attributes = c("external_gene_name", 
                              "chromosome_name",
                              "start_position",
                              "end_position"),
               mart = ensembl)

genes_clean <- genes %>%
  filter(chromosome_name %in% as.character(1:25)) %>%
  mutate(
    chromosome_name = paste0("chr", chromosome_name),
    # Add a numeric version for sorting
    chr_num = as.numeric(gsub("chr", "", chromosome_name))
  ) %>%
  arrange(chr_num)

# Create bed-like format for gene density
bed_df <- genes_clean %>%
  group_by(chromosome_name) %>%
  summarise(chr_length = max(end_position)) %>%
  ungroup()


# Calculate gene density in 2Mb windows
window_size <- 2e6
density_data <- genes_clean %>%
  group_by(chromosome_name) %>%
  do({
    # Get the range of positions
    min_pos <- min(.$start_position)
    max_pos <- max(.$end_position)
    
    # Create windows that fully span the range
    windows <- seq(floor(min_pos/window_size) * window_size,
                   ceiling(max_pos/window_size) * window_size,
                   by = window_size)
    
    # Calculate counts
    counts <- hist(.$start_position, 
                   breaks = windows, 
                   plot = FALSE)$counts
    
    data.frame(
      start = windows[-length(windows)],
      end = windows[-1],
      density = counts
    )
  }) %>%
  ungroup()

# Export gene density data
write.csv(density_data, file = "intermediate_files/20250115-gene_denisty_GRCz11.csv", row.names = FALSE)

# In case Ensembl is unresponsive use the following precomputed gene density:
#gene_density <- read.csv("/Users/ferenc.kagan/Documents/Projects/fish_transposon/intermediate_files/20250115-gene_denisty_GRCz11.csv")

gene_density <- gene_density[, -1]

# Generate the circos plot for Tuebingen
strain_data <- strain_results[["Tuebingen"]]
genome_bins <- strain_data$genome_bins
te_bins <- strain_data$te_bins %>%
  filter(!(is.na(TE))) %>%
  mutate(TE = factor(TE, levels = c("EnSpm-N49", "EnSpm-N49B", "EnSpm-N49/N49B")))

te_bins <- genome_bins %>%
  distinct(`Chromosome name`, Start, End, Bin) %>%  # Get unique chromosome bins
  expand(
    `Chromosome name`, Start, Bin, TE=unique(te_bins$TE)) %>% 
  mutate(End = Start + 2e6 - 1) %>% 
  left_join(te_bins, by = c("Chromosome name", "Start", "End", "Bin", "TE")) %>% 
  mutate(hit_count = ifelse(is.na(hit_count), 0, hit_count))

# Define xlim for Tuebingen chromosomes
chromosome_levels <- paste0("chr", 1:25)
genome_bins$`Chromosome name` <- factor(genome_bins$`Chromosome name`, levels = paste0("chr", 1:25), ordered = T)
xlim <- genome_bins %>%
  group_by(`Chromosome name`) %>%
  summarize(Start = min(Start), End = max(End)) %>%
  ungroup() %>%
  dplyr::select(Start, End) %>%
  as.matrix()

# Create circos plot for Tuebingen
pdf("output/circos_Tuebingen_with_gene_density.pdf", width = 9, height = 9)
circos.clear()

circos.par(
  start.degree = 90, 
  gap.degree = c(rep(1, length(chromosome_levels) - 1), 4),
  points.overflow.warning = FALSE)

circos.initialize(
  factors = chromosome_levels,
  xlim = cbind(xlim[, 1], xlim[, 2])
)

# Prepare colors for TEs
te_colors <- setNames(c("#FFAC94", "#FF83B1", "#7FAF75"), 
                        levels(te_bins$TE))

# Add TE hit tracks
for (te in levels(te_bins$TE)) {
  te_data <- te_bins %>% 
    filter(TE == !!te) %>%
    dplyr::select(-Bin) %>%
    distinct(.keep_all = T)
  
  circos.trackPlotRegion(
    factors = te_data$`Chromosome name`,
    y = te_data$hit_count,
    x = te_data$Start,
    track.height = 0.1,
    ylim = c(0, 4),
    panel.fun = function(region_start, y) {
      circos.barplot(y, region_start, col = te_colors[[te]], border = te_colors[[te]])
    }
  )
}

set_track_gap(cm_h(0.8))

# Add gene density track (line plot)
gene_density$`Chromosome name` <- factor(
  gene_density$chromosome_name,
  levels = chromosome_levels,
  ordered = TRUE
)

circos.trackPlotRegion(
  factors = gene_density$`Chromosome name`,
  y = gene_density$density,
  x = gene_density$start,
  track.height = 0.1,
  panel.fun = function(region_start, y) {
    circos.xaxis(
      #h = 'bottom',
      labels.facing = 'clockwise',
      labels.niceFacing = TRUE,
      major.at = c(0, 20e6, 40e6, 60e6),  # Positions in base pairs
      major.tick.length = 1,  # Adjust the tick size
      labels = c("0 Mb", "20 Mb", "40 Mb", "60 Mb"),  # Labels in Mb
      labels.cex = 0.4)  # Font size for labels
    circos.lines(region_start, y, col = "gray", lwd = 1, type = 's', area = TRUE)
  }
)

# Add labels
for (chromosome in unique(genome_bins$`Chromosome name`)) {
  mid_point <- genome_bins %>%
    filter(`Chromosome name` == chromosome) %>%
    summarize(mid = mean(c(Start, End)))
  
  circos.text(
    sector.index = chromosome,
    track.index = 1,
    x = mid_point$mid,
    y = max(te_bins$hit_count) * 1.5,  # Adjust y for placement
    labels = chromosome,
    cex = 0.8,
    facing = "bending.inside"
  )
}

text(0, 0, "Strain: Tuebingen", cex = 1.2)
legend(x = -0.23, y = -0.1,
       pch = 16, 
       legend = names(te_colors), 
       col = te_colors, 
       pt.cex = 1.5, 
       cex = 1, 
       text.col = "black", 
       horiz = F, 
       inset = c(0.1, 0.1))

# Y-axis labels
circos.yaxis(sector.index = "chr1",
               at = c(seq(1, 4, 1)), 
               track.index = 1,
               labels.cex = 0.5)
circos.yaxis(sector.index = "chr1",
               at = c(seq(1, 4, 1)), 
               track.index = 2,
               labels.cex = 0.5)
circos.yaxis(sector.index = "chr1",
               at = c(seq(1, 4, 1)), 
               track.index = 3,
               labels.cex = 0.5)

circos.yaxis(sector.index = "chr1",
               track.index = 4,
               labels.cex = 0.3)


circos.clear()
dev.off()


###########################################
### Annotate transposon insertion sites ###
###########################################

transp <- hits %>% 
  filter(strain == "Tuebingen",
         str_detect(`UCSC style name`, "chr[0-9]+$")) %>%
  rename("query" = "query acc.ver")


for (transposon in unique(transp$query)) {

    transp_temp <- transp %>%
        filter(query == transposon)

    transp_gr <- GRanges(seqnames = transp_temp$`UCSC style name`,
                         ranges = IRanges(start = transp_temp$`s. start`,
                                          end = transp_temp$`s. end`))

    # Annotate transposon insertion sites
    peaks <- annotatePeak(
        transp_gr,
        TxDb = txdb,
        level = "transcript",
        assignGenomicAnnotation = TRUE,
        annoDb = NULL,
        addFlankGeneInfo = FALSE,
        flankDistance = 5000,
        sameStrand = FALSE,
        ignoreOverlap = FALSE,
        ignoreUpstream = FALSE,
        ignoreDownstream = FALSE,
        overlap = "TSS",
        verbose = TRUE
    )

    if(transposon == "EnSpm-N49/N49B"){
      transposon <- "EnSpm-N49_N49B"
    }
    # Create and save an upset plot for annotation results
    output_path <- paste0(
        "../output/transposon_coordinate_annotations_upsetplot_", 
        transposon,
        ".pdf"
    )

    pdf(output_path, width = 10, height = 10)
    print(upsetplot(peaks, vennpie = F))
    dev.off()

    output_path <- paste0(
    "../output/transposon_coordinate_annotations_vennpie_", 
    transposon,
    ".pdf"
    )

    pdf(output_path, width = 10, height = 10)
    vennpie(peaks)
    dev.off()

}

###################################################
###                                             ###
### Tandem Site Duplication Enrichment Analysis ###
###                                             ###
###################################################
# Join BLAST output with expected hit lengths and compute match proportion
blastout <- blastout %>%
  left_join(expect, by = c("X1" = "name")) %>%
  mutate(proportion = (X8 - X7) / length) %>%  # Calculate proportion of query covered by hit
  filter(proportion >= 0.8) %>%  # Keep only high-confidence hits (>= 80% coverage)
  mutate(X2 = paste0("chr", 
                     chr_crossref$`Chromosome name`[match(.$X2, chr_crossref$`RefSeq seq accession`)]))  # Replace accession with chromosome name

# Rename BLAST output columns for readability
colnames(blastout) <- c(
  "query acc.ver", "subject acc.ver", "% identity", "alignment length",
  "mismatches", "gap opens", "q. start", "q. end",
  "s. start", "s. end", "evalue", "bit score"
)

# Create genomic ranges for Tuebingen strain hits
te_granges <- GRanges(
  seqnames = blastout$`subject acc.ver`,
  IRanges(
    start = pmin(blastout$`s. start`, blastout$`s. end`),  # Ensure start is always <= end
    end = pmax(blastout$`s. start`, blastout$`s. end`)
  )
)

# Read the reference genome (limit to first 25 chromosomes)
genome <- readDNAStringSet("../fasta/GCA_033170195.3_ASM3317019v3_genomic.fna")
genome <- genome[1:25]
names(genome) <- str_replace_all(str_extract_all(names(genome), "chromosome [0-9]+"), "chromosome ", "chr")

# Compute dinucleotide frequencies for the whole genome (background)
genome_freq <- oligonucleotideFrequency(genome, width = 2, as.prob = TRUE) %>%
  as_tibble() %>%
  pivot_longer(cols = everything(), names_to = 'Dinucleotide', values_to = 'Frequency') %>%
  mutate(Category = "Genome\nbackground")

# Compute dinucleotide frequencies in TE-flanking regions and compare with genome background
oligonucleotideFrequency(te_sequences, width = 2, as.prob = TRUE) %>%
  as_tibble() %>%
  pivot_longer(cols = everything(), names_to = 'Dinucleotide', values_to = 'Frequency') %>%
  mutate(Category = "TE\nflanking") %>%
  rbind(., genome_freq) %>%
  tidyplot(x = Category, y = Frequency, color = Category) |> 
  add_data_points_jitter(white_border = TRUE, alpha = 0.3) |>  # Add jittered points for visibility
  add_boxplot() |>  # Overlay boxplot
  add_test_pvalue(p.adjust.method = 'BH', 
                  hide.ns = TRUE, 
                  hide_info = TRUE, 
                  label = "p.adj.signif",
                  method = "wilcox_test") |>  # Wilcoxon test for pairwise comparison
  remove_x_axis_title() |>  # Clean up plot aesthetics
  remove_legend() |> 
  adjust_font(fontsize = 10) |>
  split_plot(by = Dinucleotide, height = 30, width = 50) |>  # Facet plot by dinucleotide
  save_plot("../output/dinucl_freq.pdf")  # Save plot

# Generate genomic ranges for TE flanking regions (±50 bp around TE hit)
flanking_size <- 50
te_flanks <- resize(te_granges, GenomicRanges::width(te_granges) + 2 * flanking_size, fix = "center")

# Extract TE-flanking sequences from the genome
te_sequences <- genome[te_flanks]

# Run motif discovery using MEME on TE-flanking sequences
meme_results <- run_meme(
  te_sequences,
  bin = "/Users/ferenc.kagan/Documents/Tools/meme/bin/meme",  # Path to MEME binary
  minw = 3,       # Minimum motif width
  maxw = 9,       # Maximum motif width
  nmotif = 5,     # Number of motifs to discover
  revcomp = TRUE, # Search reverse complements
  mod = "zoops",  # Zero or one occurrence per sequence
  objfun = "de"   # Objective function for discrimination
)

# View top discovered motifs (metadata)
head(meme_results)

# Plot sequence logos for each motif and save as PNG
for (i in seq_along(meme_results$motifs)) {
  motif_pwm <- convert_motifs(meme_results$motifs[[i]], class = "seqLogo-pwm")  # Convert motif to PWM format
  
  png(paste0("../output/seqlogo", i, ".png"), width = 4, height = 3, units = 'in', res = 300)
  seqLogo(motif_pwm, xaxis = FALSE, yaxis = TRUE)  # Plot logo
  dev.off()
}

# Define paths to saved motif logo images
path_images <- c(
  "/Users/ferenc.kagan/Documents/Projects/fish_transposon/output/seqlogo1.png",
  "/Users/ferenc.kagan/Documents/Projects/fish_transposon/output/seqlogo2.png",
  "/Users/ferenc.kagan/Documents/Projects/fish_transposon/output/seqlogo3.png",
  "/Users/ferenc.kagan/Documents/Projects/fish_transposon/output/seqlogo4.png",
  "/Users/ferenc.kagan/Documents/Projects/fish_transposon/output/seqlogo5.png"
)

# Create data frame summarizing motifs with consensus sequences and logos
motif_df <- data.frame(
  Name = paste0("Motif-", 1:5),
  Consensus = sapply(meme_results$motifs, function(x) x@consensus),  # Extract consensus sequence
  E_value = formatC(sapply(meme_results$motifs, function(x) x@eval), format = "e", digits = 1),  # Format E-values
  Logo = path_images,  # Paths to logos
  stringsAsFactors = FALSE
)

# Create nicely formatted table with embedded motif logos
ft <- motif_df %>%
  flextable() %>%
  compose(
    j = "Logo",
    value = as_paragraph(
      as_image(src = Logo, width = 1.5, height = 1, unit = "in")  # Embed images
    )
  ) %>%
  set_header_labels(
    Name = "Motif Name",
    Consensus = "Consensus Sequence",
    E_value = "E-value",
    Logo = "Sequence Logo"
  ) %>%
  fontsize(size = 11, part = "all") %>%
  align(align = "center", part = "all") %>%
  autofit()

# Display flextable in interactive R session or export using `save_as_docx(ft)` or similar
ft

