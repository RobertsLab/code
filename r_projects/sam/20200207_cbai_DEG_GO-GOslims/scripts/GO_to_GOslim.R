library(GSEABase)
library(tidyverse)

#########################################################################
# Scipt to map C.bairdi DEG enriched GO terms to GOslims
# and identify the GO terms contributing to each GOslim
#########################################################################



### Set false discovery rate (FDR) filter, if desired
fdr <- as.character("1.0")

### Download files
download.file(url = "https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/infected-vs-uninfected/edgeR.2317.dir/salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.infected-UP.subset.GOseq.enriched.flattened",
              destfile = "./data/infected-vs-uninfected/salmon.gene.counts.matrix.infected_vs_uninfected.edgeR.DE_results.P0.05_C1.infected-UP.subset.GOseq.enriched.flattened")
download.file(url = "https://gannet.fish.washington.edu/Atumefaciens/20200207_cbai_DEG/D12-vs-D26/edgeR.21229.dir/salmon.gene.counts.matrix.D12_vs_D26.edgeR.DE_results.P0.05_C1.D26-UP.subset.GOseq.enriched.flattened", 
              destfile = "./data/D12-vs-D26/salmon.gene.counts.matrix.D12_vs_D26.edgeR.DE_results.P0.05_C1.D26-UP.subset.GOseq.enriched.flattened")


### Create list of files
goseq_files <- list.files(path = "./data",
                         pattern = "\\.GOseq.[de]",
                         recursive = TRUE,
                         full.names = TRUE)

### Set output filename suffix
output_suffix=("GOslims.csv")

### Strip path from goseq files
goseq_filename <- basename(goseq_files)

### Vector of GOslim ontologies (e.g. Biological Process = BP, Molecular Function = MF, Cellular Component = CC)
ontologies <- c("BP", "CC", "MF")

for (slim_ontology in ontologies) {

  ### Set GOOFFSPRING database, based on ontology group set above
  go_offspring <- paste("GO", slim_ontology, "OFFSPRING", sep = "")
  
  for (item in goseq_files) {
    
    ## Get max number of fields
    # Needed to handle reading in file with different number of columns in each row
    # as there may be multiple 
    max_fields <- max(count.fields(item, sep = "\t"), na.rm = TRUE)
    
    ## Read in tab-delimited GOseq file
    # Use "max_fields" to populate all columns with a sequentially numbered header
    go_seqs <- read.delim(item,
                          col.names = paste0("V",seq_len(max_fields)))
    
    ## Filter enriched GOterms with false discovery rate
    goseqs_fdr <- filter(go_seqs, V8 <= as.numeric(fdr))
    
    ## Grab just the individual GO terms from the "category" column)
    goterms <- as.character(goseqs_fdr$V1)
    
    ### Use GSEA to map GO terms to GOslims
    
    ## Store goterms as GSEA object
    myCollection <- GOCollection(goterms)
    
    ## Use generic GOslim file to create a GOslim collection
    
    # I downloaded goslim_generic.obo from http://geneontology.org/docs/go-subset-guide/
    # then i moved it to the R library for GSEABase in the extdata folder
    # in addition to using the command here - I think they're both required.
    slim <- getOBOCollection("./data/goslim_generic.obo")
    
    ## Map GO terms to GOslims and select Biological Processes group
    slimsdf <- goSlim(myCollection, slim, slim_ontology)
    
    ## Need to know the 'offspring' of each term in the ontology, and this is given by the data in:
    # GO.db::getFromNamespace(go_offspring, "GO.db")
    
    ## Create function to parse out GO terms assigned to each GOslim
    ## Courtesy Bioconductor Support: https://support.bioconductor.org/p/128407/
    mappedIds <-
      function(df, collection, OFFSPRING)
      {
        map <- as.list(OFFSPRING[rownames(df)])
        mapped <- lapply(map, intersect, ids(collection))
        df[["go_terms"]] <- vapply(unname(mapped), paste, collapse = ";", character(1L))
        df
      }
    
    ## Run the function
    slimsdf <- mappedIds(slimsdf, myCollection, getFromNamespace(go_offspring, "GO.db"))
    
    ## Provide column name for first column
    slimsdf <- cbind(GOslim = rownames(slimsdf), slimsdf)
    rownames(slimsdf) <- NULL
    
    ### Prep output file naming structure
    
    ## Split filenames on periods
    ## Creates a list
    split_filename <- strsplit(item, ".", fixed = TRUE)
    
    ## Split filename on directories
    ## Creates a list
    split_dirs <- strsplit(item, "/", fixed = TRUE)
    
    ## Slice split_filename list from position 9 to the last position of the list
    ## Paste these together using a period
    goseq_filename_split <-paste(split_filename[[1]][9:lengths(split_filename)], collapse = ".")
    
    ## Slice split_dirs list at position
    
    ## Paste elements together to form output filename
    fdr_file_out <- paste("FDR", fdr, sep = "_")
    outfilename <- paste(goseq_filename_split, fdr_file_out, slim_ontology ,output_suffix, collapse = ".", sep = ".")
    
    ## Set output file destination and name
    ## Adds proper subdirectory from split_dirs list
    outfile_dest <- file.path("./analyses", split_dirs[[1]][3], outfilename)
    
    ## Write output file
    write.csv(slimsdf, file = outfile_dest, quote = FALSE, row.names = FALSE)
    
    
  }
}


