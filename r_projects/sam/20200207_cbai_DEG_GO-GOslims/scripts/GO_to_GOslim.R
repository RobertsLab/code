library(GSEABase)
library(tidyverse)

download.file(url = "https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/infected-vs-uninfected/edgeR.2317.dir/salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.infected-UP.subset.GOseq.enriched",
              destfile = "./data/salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.infected-UP.subset.GOseq.enriched")
download.file(url = "https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/infected-vs-uninfected/edgeR.2317.dir/salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.infected-UP.subset.GOseq.depleted",
              destfile = "./data/salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.infected-UP.subset.GOseq.depleted")
download.file(url = "https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/infected-vs-uninfected/edgeR.2317.dir/salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.uninfected-UP.subset.GOseq.enriched",
              destfile = "./data/salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.uninfected-UP.subset.GOseq.enriched")
download.file(url = "https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/infected-vs-uninfected/edgeR.2317.dir/salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.uninfected-UP.subset.GOseq.depleted",
              destfile = "./data/salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.uninfected-UP.subset.GOseq.depleted")

goseq_files <- list.files(path = "./data",
                         pattern = "\\.GOseq.[de]",
                         full.names = TRUE)
  
output_suffix=("GOslims.csv")



# Set full paths for goseq files
goseq_filename <- basename(goseq_files)
  
for (item in goseq_files) {
  
  print(item)
  
  max_fields <- max(count.fields(item, sep = "\t"))
  
  # Read in tab-delimited GOseq file
  go_seqs <- read.table(item,
                        sep = "\t",
                        header = TRUE,
                        col.names = paste0("V",seq_len(max_fields)),
                        fill = TRUE)
  
  # Grab just the individual GO terms (i.e. "category" column)
  goterms <- as.character(go_seqs$category)
  
  #goslims with GSEA
  myCollection <- GOCollection(goterms)
  
  #I downloaded goslim_generic.obo from http://geneontology.org/docs/go-subset-guide/
  #then i moved it to the R library for GSEABase in the extdata folder
  slim <- getOBOCollection("~/Downloads/goslim_generic.obo")
  
  # Select Biological Processes groups
  slims <- goSlim(myCollection, slim, "BP")
  
  
  ### Prep output file naming structure
  
  ## Split filenames on periods
  ## Creates a list
  split_filename <- strsplit(item, ".", fixed = TRUE)

  ## Slice split_filename list from position 9 to the last position of the list
  ## Paste these together using a period
  goseq_filename_split <-paste(split_filename[[1]][9:lengths(split_filename)], collapse = ".")
  
  ## Paste elements together to form output filename
  outfilename <- paste(goseq_filename_split, output_suffix, collapse = ".", sep = ".")
  
  # Set output file destination and name
  outfile_dest <- file.path("./analyses/", outfilename)
  
  ## Write output file
  write.csv(slims, file = outfile_dest, quote = FALSE)
  
  
}



