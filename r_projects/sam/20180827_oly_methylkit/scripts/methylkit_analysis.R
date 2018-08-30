######
# This script uses the MethylKit package to examine methylation
# differences between a single oyster population raised in two different locations
# in Puget Sound, WA, USA - Fidalgo Bay & Oyster Bay.
# 
# The analysis is set to analyze CpG having minimum 3x coverage.
# This can be changed in the processBismarkAln function below
######


library(tidyverse)
library(methylKit)

# Download deduplicated, sorted, BAM files
download.file("http://owl.fish.washington.edu/Athaliana/20180709_oly_methylseq/1_ATCACG_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/1_ATCACG_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180709_oly_methylseq/2_CGATGT_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/2_CGATGT_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180709_oly_methylseq/3_TTAGGC_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/3_TTAGGC_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180709_oly_methylseq/4_TGACCA_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/4_TGACCA_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180709_oly_methylseq/5_ACAGTG_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/5_ACAGTG_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180709_oly_methylseq/6_GCCAAT_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/6_GCCAAT_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180709_oly_methylseq/7_CAGATC_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/7_CAGATC_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180709_oly_methylseq/8_ACTTGA_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/8_ACTTGA_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")

# Store bam files in list
bam_files_list <- as.list(list.files(path = "./data",
                             pattern = "\\.bam$",
                             full.names = TRUE))

# List of sample IDs
## 1 = Fidalgo Bay outplant
## 2 = Oyster Bay outplant
## NF = Fidalgo Bay broodstock origination
sample_ids_list <- list("1NF11", "1NF15", "1NF16", "1NF17",
                    "2NF5", "2NF6", "2NF7", "2NF8")

# Set number of 1NF or 2NF treatments: 0 or 1
treatmentSpecification <- c(rep(0, times = 4), rep(1, times = 4))

# Get methylation stats for CpGs with at least 3x coverage
meth_stats <- processBismarkAln(location = bam_files_list,
                                    sample.id = sample_ids_list,
                                    assembly = "Olurida_v080.fa ",
                                    read.context = "CpG",
                                    mincov = 3,
                                    treatment = treatmentSpecification)
# File count
nFiles <- length(bam_files_list)


# Generate and save histograms showing Percent CpG Methylation

for(i in 1:nFiles) {
  cpg_methylation_percent_path <- file.path("./analyses", paste("cpg_methylation_percent_", sample_ids_list[[i]], ".png", sep = "")) #Specify save destination and filename
  png(cpg_methylation_percent_path, height = 1000, width = 1000) #Save file with designated name
  getMethylationStats(meth_stats[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device
}

# Generate and save histograms showing CpG Methylation Coverage
for(i in 1:nFiles) { #For each data file
  cpg_coverage_path <- file.path("./analyses/", paste("cpg_coverage_", sample_ids_list[[i]], ".png", sep = "")) #Specify save destination and filename
  png(cpg_coverage_path, height = 1000, width = 1000) #Save file with designated name
  getCoverageStats(meth_stats[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device
}

# Combine all samples to evaluate bases covered in all samples
methylation_information <- unite(meth_stats)

# Clustering dendrogram
dendrogram_path <- file.path("./analyses", paste("clustering_dendrogram", ".png", sep = "")) #Specify save destination and filename
png(dendrogram_path, height = 1000, width = 1000) #Save file with designated name
clusterSamples(methylation_information, dist="correlation", method="ward", plot=TRUE)
dev.off()

# Run a PCA analysis on percent methylation for all samples
pca_path <- file.path("./analyses", paste("pca", ".png", sep = "")) #Specify save destination and filename
png(pca_path, height = 1000, width = 1000) #Save file with designated name
PCASamples(methylation_information)
dev.off()

#Run the PCA analysis and plot variances against PC number in a screeplot
scree_path <- file.path("./analyses/", paste("pca_scree", ".png", sep = "")) #Specify save destination and filename
png(scree_path, height = 1000, width = 1000) #Save file with designated name
PCASamples(methylation_information, screeplot = TRUE)

#Calculate differential methylation statistics based on treatment indication from processBismarkAln
differentialMethylationStats <- calculateDiffMeth(methylation_information)

#Identify loci that are at least 25% different. Q-value is the FDR used for p-value corrections.
diffMethStats25 <- getMethylDiff(differentialMethylationStats)

#Convert to bedgraph
bedgraph(diffMethStats25, file.name = "./analyses/oly_fb_ob_dml.bed", col.name = "meth.diff")

