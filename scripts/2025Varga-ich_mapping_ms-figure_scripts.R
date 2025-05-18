# Scripts for generating some figures for the Varga et al (2025) manuscript
# on the mapping of the zebrafish ichabod (ich) mutation

library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggpubr)


# Dotplots showing the alighnments between ONT reads and chromosome 19 
# Figure 1A

library(pafr)
library(tidyverse)

ich_ch19_al <- read_paf("all_ich_cntg_vs_TU_chr19_Fish6.paf")

ggplot(ich_ch19_al, aes(alen)) + 
  geom_histogram(colour="black", fill="steelblue", bins=20) + 
  theme_bw(base_size=16) + 
  ggtitle("Distribution of alignment lengths") +
  scale_x_log10("Alignment-length")

long_ich_ch19_al <- subset(ich_ch19_al, alen > 2e5)

dotplot(long_ich_ch19_al, label_seqs=FALSE) + 
  theme_bw() + coord_flip()

ggsave("20250104_TU_chr19_Fish6_vs_all_ich_contigs_dotplot.pdf", 
       units="cm", width=15, height=15)

# Dotplots showing the alighnments between ONT conti tig00035647 and the telomeric end of chromosome 19 
# Figure 1B (upper part)

plot_synteny(long_ich_ch19_al, 
             q_chrom="tig00035647", 
             t_chrom="chr19", 
             centre=FALSE
) +
  theme_bw()

plot_coverage(long_ich_ch19_al)   

ich_ctnnb2_al <- read_paf("ich_ctnnb2_vs_TU_ch19_Fish6_end.paf")

ich_ctnnb2_al  <- subset(ich_ctnnb2_al, mapq >31)
ich_ctnnb2_al  <- subset(ich_ctnnb2_al, alen >1e4)

dotplot(ich_ctnnb2_al, label_seqs=FALSE) + 
  theme_bw() + coord_flip()

ggsave("20250104_TU_chr19_Fish6_end_vs_ctnnb2_ich_contigs_dotplot.pdf", 
       units="cm", width=10, height=10)


# Plotting the genes in the relevant region of chromosome 19
# Figure 1B (lower part)

library(rtracklayer)
library(GenomicRanges)
library(Gviz)

# Import GTF file
gtf <- import("NHGRI_Fish6.gtf", format="GTF")

# Create tracks
gtrack <- GenomeAxisTrack()
grtrack <- GeneRegionTrack("NHGRI_Fish6.gtf",
                           chromosome="chr19",
                           #geneSymbol = TRUE,
                           #transcriptAnnotation = "gene_name",
                           #name = "Gene" 
)


antrack <- AnnotationTrack("NHGRI_Fish6.gtf",
                           chromosome="chr19",
                           name="Genes",
                           shape="box",
                           fill="red")

# Plot
pdf("20250104_TU_chr19_Fish6_end_genes.pdf")
plotTracks(list(gtrack, antrack),
           from=49487817,
           to=51759963,
           chromosome="chr19",
           collapseTranscripts = "meta", 
           shape = "arrow", 
           transcriptAnnotation = "symbol",
           #fill="darkred",
           #color = "black",
           sizes=c(1,1))

dev.off()

ggsave("20250104_TU_chr19_Fish6_end_genes.pdf", 
       units="cm", width=10, height=5)


# Reanalysis of SLAM-Seq data-sets provided in the Supplementary Information of Bhat et al., 2023  
# Figures 2A and B

gene_expr_data <- read.csv("SLAMseq_reanalysis/bhat2023-SLAMseq_expr_data.csv")

conversion_data <- read.csv("SLAMseq_reanalysis/bhat2023-SLAMseq-TC_conversion.csv")

ctnnb_expr_data <- gene_expr_data %>% select(-classification) %>% select(-sub.classification) %>%
  filter(Gene.ID %in% c("ctnnb1", "ctnnb2")) %>%
  pivot_longer(cols = (-"Gene.ID"), names_to=c("blank", "rep", "time_point"), 
               names_sep="_", 
               values_to="rpm") %>%
  select(-blank)

expression <- ggplot(ctnnb_expr_data, aes(x=time_point, y=rpm, fill=Gene.ID)) +
  geom_boxplot(position=position_dodge(1)) + theme_classic() +
  scale_fill_brewer(palette = "Paired")

expression_bar <- ggplot(ctnnb_expr_data, aes(x=time_point, y=rpm, fill=Gene.ID)) +
  geom_bar(stat = "identity", position="fill") + theme_classic() +
  scale_fill_brewer(palette = "Paired")

ctnnb_conversion_data <- conversion_data %>% select(-MZ_cluster) %>%
  filter(gene %in% c("ctnnb1", "ctnnb2")) %>%
  pivot_longer(cols = (-"gene"), names_to=c("blank1", "blank2", "time_point"), 
               names_sep="_", 
               values_to="fractionTC") %>%
  select(-blank1, -blank2) 

conversion <-ggplot(ctnnb_conversion_data, aes(x=time_point, y=fractionTC, color=gene)) +
  geom_point() + theme_classic() +
  scale_color_brewer(palette = "Paired")

ggarrange(expression_bar, conversion, common.legend = TRUE)

ggsave("20240116-ctnnb1&2_SLAMseq.pdf", units = "cm", width=20, height = 10)

# Boxplots showing the results of relative fluorescence measurements for 3'UTR stability 
# Figure 2G

UTR3_stability_data <- read.csv("ich_3UTR_stability_results-20250110.csv")

UTR3_stability_data <- UTR3_stability_data %>% group_by(gene) %>%
  pivot_wider(names_from = "marker", values_from = "mean.exp") %>%
  mutate(ratio = gfp/mch) %>% filter(gene!="h2afx1") %>%
  filter(sample_no != "ctnnb2_ich_20") %>%
  mutate(gene = factor(gene, 
                       levels= c("ccna1","ctnnb2", "cntnnb2-ich", 
                                 "tuba8l")))

my_comparisons <- list( c("ctnnb2", "cntnnb2-ich"), c("ccna1", "tuba8l"), 
                        c("ctnnb2", "ccna1"), c("cntnnb2-ich", "tuba8l"))

ggplot(UTR3_stability_data, aes(x=gene, y=ratio, fill=gene)) +
  geom_boxplot(position=position_dodge(1)) + theme_classic() +
  scale_fill_manual(values = c("#E7EA8E", "#1F78B4", "#7FAF75", "#BA6BA9")) + 
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value


ggsave("20250110-ich_3UTR_GFP_stability.pdf", units = "cm", width=10, height=10)

# Protein structures for the two zebrafish beta catenin proteins drawn with the
# drawProteins Bioconductor package (DOI: 10.18129/B9.bioc.drawProteins)
# Supplementary Figure 2A

if (!require(devtools)) {
  install.packages('devtools')
}
dev_mode(on=TRUE)
devtools::install_github('brennanpincardiff/drawProteins')

library(drawProteins)
library(ggplot2)
library(magrittr)

# accession numbers of the two beta-catenin proteins
"Q8JID2 A0A5H1ZRJ2" %>%
  drawProteins::get_features() %>%
  drawProteins::feature_to_dataframe() ->
  prot_data

p <- draw_canvas(prot_data)
p <- draw_chains(p, prot_data, label_size = 2)
p <- draw_domains(p, prot_data)
p <- draw_repeat(p, prot_data, fill = "darkgrey")
p <- draw_regions(p, prot_data)

#p <- draw_phospho(p, prot_data, size = 8)

# background and y-axis
p <- p + theme_bw(base_size = 20) +  # white background and change text size
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  theme(#axis.ticks = element_blank(),
    axis.text.y = element_blank()) +
  theme(panel.border = element_blank())

# add titles
p <- p + labs(x = "Amino acid number",         # label x-axis
              y = "",  # label y-axis
              title = "Schematic of zebrafish beta-catenin proteins",
              #subtitle = "circles = phosphorylation sites\nRHD = Rel Homology Domain\nsource:Uniprot"
) +
  scale_x_continuous(name="Amino acid number", limits=c(0, 800))


# move legend to top
p <- p + theme(legend.position="top") + labs(fill="")

p

ggsave("20250106-beta_catenin_domain_structures.pdf", units="cm", width=20, height=12)


# Supplementary Figure 3C was generated using the raw data provided as Supplementary Information
# for the Chang et al 2022 Genome Research paper: https://genome.cshlp.org/content/32/7/1408/suppl/DC1

telescope_N49_annot <- read_table("../Chang2022_reanalysis/TE_annotations_DEclusters.gtf", 
                                  skip=3, col_names = FALSE) %>%
  filter(str_detect(X10, "^\"EnSpm-N49"))

write_tsv(telescope_N49_annot, "../Chang2022_reanalysis/EnSpm-N49_like_clusters.gtf", 
          col_names = FALSE, quote= "none")  

telescope_N49_annot <- read_table("../Chang2022_reanalysis/EnSpm-N49_like_clusters.gtf", 
                                  skip=3, col_names = FALSE) %>% 
  filter(!str_detect(X20, "^\"Not_DE"))


telescope_N49_counts <- read_table("../Chang2022_reanalysis/Telescop_TEloci_counts/Telescope_TE_counts_normalized.tab") %>%
  filter(str_detect(TEs, "^EnSpm-N49")) %>% select(!"TE2s") %>%
  gather(key = "Stage_Rep", value = "Expression", -TEs) %>%
  mutate(Stage = str_extract(Stage_Rep, "^[^_]+.*?(?=_[rR]ep_)")) %>%
  group_by(TEs, Stage) %>%
  summarize(Mean_Expression = mean(Expression))

telescope_N49_counts$Stage <- factor(telescope_N49_counts$Stage, 
                                     levels=c("1_cell", "2_cell", "128_cell", "1k_cell", "dome",
                                              "50pc_epiboly", "shield","75pc_epiboly","1_4_somites",
                                              "14_19_somites", "20_25_somites", "prim_5",
                                              "prim_15", "prim_25", "long_pec", "protruding_mouth",
                                              "day_4", "day_5"))

temp_df <- telescope_N49_counts %>% 
  mutate(TE_fam = str_extract(TEs, "^[^_]+.*?(?=_DR_)"))

# Create heatmap
ggplot(telescope_N49_counts, aes(x = Stage, y = TEs, fill = Mean_Expression)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "red") + 
  labs(x = "Stage", y = "TEs", fill = "Mean Expression") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(temp_df %>% group_by(TE_fam), aes(x = Stage, y = TEs, fill = Mean_Expression)) + 
  #facet_wrap(vars(TE_fam), nrow=2, scales = "free") +
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "red") + 
  labs(x = "Stage", y = "TEs", fill = "Mean Expression") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("20250101-EnSpm-N49_49B-heatmap.pdf", units="cm", width=20, height=20)