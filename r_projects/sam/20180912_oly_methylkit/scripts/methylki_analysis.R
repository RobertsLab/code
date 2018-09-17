######
# This script uses the MethylKit package to examine methylation
# differences between a single oyster population raised in two different locations
# in Puget Sound, WA, USA - Fidalgo Bay & Oyster Bay.
#
######


library(tidyverse)
library(methylKit)

# Saves R and R Studio version info to file.
write(capture.output(sessionInfo()), file = "version_info.txt")
write(capture.output(RStudio.Version()), file = "version_info.txt", append = TRUE)


download.file("http://owl.fish.washington.edu/Athaliana/20180912_oly_WGBSseq_bismark/1_ATCACG_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/1_ATCACG_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180912_oly_WGBSseq_bismark/2_CGATGT_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/2_CGATGT_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180912_oly_WGBSseq_bismark/3_TTAGGC_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/3_TTAGGC_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180912_oly_WGBSseq_bismark/4_TGACCA_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/4_TGACCA_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180912_oly_WGBSseq_bismark/5_ACAGTG_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/5_ACAGTG_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180912_oly_WGBSseq_bismark/6_GCCAAT_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/6_GCCAAT_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180912_oly_WGBSseq_bismark/7_CAGATC_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/7_CAGATC_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")
download.file("http://owl.fish.washington.edu/Athaliana/20180912_oly_WGBSseq_bismark/8_ACTTGA_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam", destfile = "./data/8_ACTTGA_L001_R1_001_trimmed_bismark_bt2.deduplicated.sorted.bam")

# List of sample IDs
## 1 = Fidalgo Bay outplant
## 2 = Oyster Bay outplant
## NF = Fidalgo Bay broodstock origination
fidalgo_bay_ids_list <- list("1NF11", "1NF15", "1NF16", "1NF17")

oyster_bay_ids_list <- list("2NF5", "2NF6", "2NF7", "2NF8")

# Assign treatment binaries

fidalgo_bay_treatment = 0

oyster_bay_treatment = 1

# Set number of 1NF or 2NF treatments.
# 0 or 1 is assigned to each sample ID corresponding to a particular treatment.
treatmentSpecification <- c(
  rep(fidalgo_bay_treatment, times = length(fidalgo_bay_ids_list)), 
  rep(oyster_bay_treatment, times = length(oyster_bay_ids_list)))
