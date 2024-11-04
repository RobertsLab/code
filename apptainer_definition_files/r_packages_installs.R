# Update base packages
update.packages(ask = FALSE)

# Install BioConductor package manager
if (!requireNamespace("BiocManager", quietly = TRUE, ask = FALSE))
install.packages("BiocManager")
BiocManager::install(version = "3.19")

# Install tidyverse
install.packages("tidyverse")

# Install matrixStats 0.61.0 (needed for DESeq2)
install.packages("https://cran.rstudio.com/src/contrib/matrixStats_0.61.0.tar.gz", repos=NULL, type="source")

# Install remotes package (allows for package installs from GitHub)
BiocManager::install("remotes")

# Install GSEABase (a dependency for numerous gene ontology/enrichment analysis)
BiocManager::install("Bioconductor/GSEABase")

# Install qvalue package
BiocManager::install("qvalue")

# Install GO.db (annotation maps for Gene Ontology data)
BiocManager::install("GO.db")

# Install MatrixGenerics (needed for DESeq2)
BiocManager::install("MatrixGenerics")

# Install DESeq2
BiocManager::install("DESeq2")
