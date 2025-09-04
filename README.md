# SNP-Landscape-in-Head-and-Neck-Squamous-Cell-Cancer
Complete workflow for analyzing and visualizing Single Nucleotide Polymorphisms (SNPs) from Whole Exome Sequencing (WES) data of Head and Neck Squamous Cell Carcinoma (HNSCC), based on open-source data (SRA ID: SRR32633603). Using R and Bioconductor packages like VariantAnnotation and Genomic
# Load libraries
# --------------------------
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyr)  # for data reshaping

# --------------------------
# Simulate example SNP dataset
# Replace with your real data!
# --------------------------
set.seed(123)
n_snps <- 1000
n_samples <- 50  # number of samples for genotype matrix simulation

snp_df <- data.frame(
  REF = sample(c("A", "C", "G", "T"), n_snps, replace = TRUE),
  ALT = sample(c("A", "C", "G", "T"), n_snps, replace = TRUE),
  QUAL = runif(n_snps, 10, 500),
  CHR = sample(as.character(1:22), n_snps, replace = TRUE),
  POS = sample(1:5e6, n_snps, replace = TRUE)
) %>%
  filter(REF != ALT) %>%
  mutate(
    mutation = paste0(REF, ">", ALT),
    bin = floor(POS / 1e6),
    type = ifelse(
      (REF == "A" & ALT == "G") | (REF == "G" & ALT == "A") |
        (REF == "C" & ALT == "T") | (REF == "T" & ALT == "C"),
      "Transition", "Transversion"
    )
  )

# --------------------------
# Simulate genotype matrix (0,1,2 for genotypes) for PCA
# SNPs in rows, samples in columns
# --------------------------
genotypes <- matrix(
  sample(0:2, n_snps * n_samples, replace = TRUE, prob = c(0.7, 0.2, 0.1)),
  nrow = n_snps, ncol = n_samples
)
colnames(genotypes) <- paste0("Sample", 1:n_samples)
rownames(genotypes) <- paste0("SNP", 1:n_snps)

# --------------------------
# Perform PCA on genotype matrix
# --------------------------
pca_res <- prcomp(t(genotypes), scale. = TRUE)  # transpose: samples are rows

# Create PCA dataframe for plotting
pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2]
)

# --------------------------
# Plot PCA
# --------------------------
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(color = "darkblue", size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Simulated Genotype Data",
       x = "PC1", y = "PC2")

# --------------------------
# Previous plots
# --------------------------

# 1. SNP Density Heatmap
bin_counts <- snp_df %>%
  group_by(CHR, bin) %>%
  summarise(SNP_count = n(), .groups = "drop")

ggplot(bin_counts, aes(x = bin, y = CHR, fill = SNP_count)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma") +
  labs(title = "SNP Density Heatmap (per 1Mb Region)",
       x = "Genomic Bin (Mb)", y = "Chromosome", fill = "SNP Count") +
  theme_minimal()

# 2. Pie Chart of Mutation Types
mutation_counts <- snp_df %>%
  count(mutation, name = "count")

ggplot(mutation_counts, aes(x = "", y = count, fill = mutation)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Mutation Type Proportions") +
  scale_fill_brewer(palette = "Set3")

# 3. Manhattan-style SNP Quality Plot
snp_df <- snp_df %>%
  mutate(
    CHR = factor(CHR, levels = sort(unique(CHR))),
    CHR_num = as.numeric(CHR),
    global_pos = POS + (CHR_num - 1) * 1e7
  )

ggplot(snp_df, aes(x = global_pos, y = QUAL, color = CHR)) +
  geom_point(size = 0.8, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Genome-wide SNP Quality Plot",
       x = "Pseudo-Genomic Position", y = "Quality Score") +
  theme(legend.position = "none")

# 4. Boxplot: SNP Quality by Mutation
ggplot(snp_df, aes(x = mutation, y = QUAL, fill = mutation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  theme_minimal() +
  labs(title = "SNP Quality by Mutation Type",
       x = "Mutation Type", y = "Quality Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set2")

# 5. Transition vs Transversion Plot
type_counts <- snp_df %>%
  count(type)

ggplot(type_counts, aes(x = type, y = n, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Transition" = "blue", "Transversion" = "red")) +
  theme_minimal() +
  labs(title = "Transitions vs Transversions",
       x = "Type", y = "Count") +
  theme(legend.position = "none")

# 6. Histogram of Quality Scores
ggplot(snp_df, aes(x = QUAL)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of SNP Quality Scores",
       x = "QUAL Score", y = "Frequency")

# 7. SNP Count per Chromosome
ggplot(snp_df %>% count(CHR), aes(x = CHR, y = n)) +
  geom_bar(stat = "identity", fill = "green") +
  theme_minimal() +
  labs(title = "SNP Counts per Chromosome",
       x = "Chromosome", y = "SNP Count")

<img width="2100" height="1200" alt="Variant_counts_per_chromosome" src="https://github.com/user-attachments/assets/b8a0da9a-0dd7-486a-8fe3-991722895a97" />
<img width="704" height="354" alt="SNPs per chr" src="https://github.com/user-attachments/assets/e7b6ca5a-f8c1-4bf4-bec3-63ee7ed41dfa" />
<img width="704" height="354" alt="SNP_box plot" src="https://github.com/user-attachments/assets/1b735308-b2b6-443e-a92f-ac240b806c3e" />
<img width="704" height="354" alt="SNP (TransTv)" src="https://github.com/user-attachments/assets/cfa506ea-48e9-4959-b210-49d4d8bd1a9a" />
<img width="704" height="354" alt="PCA plot_SNP_HN Cancer" src="https://github.com/user-attachments/assets/0d0f3b93-5cb8-4f0b-a2aa-9e2c010a68fd" />
<img width="704" height="354" alt="GWAS plot_HN Cancer" src="https://github.com/user-attachments/assets/e2922ffa-4114-4a97-bab7-2d7ce7be6ccf" />
<img width="704" height="354" alt="density Heatmap plot_HN Cancer" src="https://github.com/user-attachments/assets/b5dcba78-2dcb-4b2d-93e6-1e8828bc18ec" />

