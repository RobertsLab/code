01-GO-to-GOslim
================
Sam White
2024-11-06

- [1 INTRO](#1-intro)
  - [1.1 CODE](#11-code)
    - [1.1.1 Variables](#111-variables)
    - [1.1.2 Set GSEAbase location and download
      `goslim_generic.obo`](#112-set-gseabase-location-and-download-goslim_genericobo)
    - [1.1.3 Read in gene/GO file](#113-read-in-genego-file)
    - [1.1.4 Remove rows with NA, remove whitespace in GO IDs column and
      keep just gene/GO IDs
      columns](#114-remove-rows-with-na-remove-whitespace-in-go-ids-column-and-keep-just-genego-ids-columns)
    - [1.1.5 “Flatten” gene and GO IDs](#115-flatten-gene-and-go-ids)
    - [1.1.6 Group](#116-group)
    - [1.1.7 Map GO IDs to GOslims](#117-map-go-ids-to-goslims)
    - [1.1.8 Extract GOslims from OBO](#118-extract-goslims-from-obo)
    - [1.1.9 Get Biological Process (BP) GOslims associated with
      provided GO
      IDs.](#119-get-biological-process-bp-goslims-associated-with-provided-go-ids)
    - [1.1.10 Map GO to GOslims](#1110-map-go-to-goslims)
    - [1.1.11 Flatten GOslims](#1111-flatten-goslims)
  - [1.2 Write slimdf to file](#12-write-slimdf-to-file)
    - [1.2.1 Sort and select slmidf](#121-sort-and-select-slmidf)
    - [1.2.2 Write formatted slim.count.df to
      file](#122-write-formatted-slimcountdf-to-file)
  - [1.3 Count unique Biological Process GO
    IDs](#13-count-unique-biological-process-go-ids)
- [2 REFERENCES](#2-references)

# 1 INTRO

This notebook is setup to [take GO IDs and categorize them into their
corresponding
GOslims](https://github.com/RobertsLab/resources/issues/2012) (GitHub
Issue). Specifically, it will use the following tab-delimited input
file:

- [der_go.tsv](https://gannet.fish.washington.edu/seashell/bu-github/project-pycno-multispecies-2023/output/06-Go-3-species/der_go.tsv)

[Desired output
format](https://github.com/RobertsLab/resources/issues/1989#issuecomment-2448757444)
(GitHub Issue comment) is:

    Term                                                    Count
    anatomical structure development                        1137
    signaling                                               600
    cell differentiation                                    533
    immune system process                                   331
    lipid metabolic process                                 212

This was performed using R, with the following packages:

- GSEABase (Martin Morgan 2017)
- GO.db (Carlson 2017)
- tidyverse (Wickham et al. 2019)

It also relied on gene ontology (GO) information from the [Gene Ontology
Consortium](https://geneontology.org/) (Ashburner et al. 2000;
Aleksander et al. 2023).

## 1.1 CODE

``` r
library(GSEABase)
```

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: annotate

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: XML

    ## Loading required package: graph

    ## 
    ## Attaching package: 'graph'

    ## The following object is masked from 'package:XML':
    ## 
    ##     addNode

``` r
library(GO.db)
```

    ## 

``` r
library(knitr)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::%within%() masks IRanges::%within%()
    ## ✖ stringr::boundary()   masks graph::boundary()
    ## ✖ dplyr::collapse()     masks IRanges::collapse()
    ## ✖ dplyr::combine()      masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ dplyr::desc()         masks IRanges::desc()
    ## ✖ tidyr::expand()       masks S4Vectors::expand()
    ## ✖ dplyr::filter()       masks stats::filter()
    ## ✖ dplyr::first()        masks S4Vectors::first()
    ## ✖ dplyr::lag()          masks stats::lag()
    ## ✖ ggplot2::Position()   masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()       masks IRanges::reduce()
    ## ✖ dplyr::rename()       masks S4Vectors::rename()
    ## ✖ lubridate::second()   masks S4Vectors::second()
    ## ✖ lubridate::second<-() masks S4Vectors::second<-()
    ## ✖ dplyr::select()       masks AnnotationDbi::select()
    ## ✖ dplyr::slice()        masks IRanges::slice()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
knitr::opts_chunk$set(
  echo = TRUE,         # Display code chunks
  eval = FALSE,        # Evaluate code chunks
  warning = FALSE,     # Hide warnings
  message = FALSE,     # Hide messages
  comment = ""         # Prevents appending '##' to beginning of lines in code output
)
```

### 1.1.1 Variables

``` r
# Column names corresponding to gene name/ID and GO IDs
GO.ID.column <- "Gene.Ontology.IDs"
gene.ID.column <- "X1"

# Relative path or URL to input file
input.file <- "https://gannet.fish.washington.edu/seashell/bu-github/project-pycno-multispecies-2023/output/06-Go-3-species/der_go.tsv"


##### Official GO info - no need to change #####
goslims_obo <- "goslim_generic.obo"
goslims_url <- "http://current.geneontology.org/ontology/subsets/goslim_generic.obo"
```

### 1.1.2 Set GSEAbase location and download `goslim_generic.obo`

``` r
# Find GSEAbase installation location
gseabase_location <- find.package("GSEABase")

# Load path to GOslim OBO file
goslim_obo_destintation <- file.path(gseabase_location, "extdata", goslims_obo, fsep = "/")

# Download the GOslim OBO file
download.file(url = goslims_url, destfile = goslim_obo_destintation)

# Loads package files
gseabase_files <- system.file("extdata", goslims_obo, package="GSEABase")
```

### 1.1.3 Read in gene/GO file

``` r
full.gene.df <- read.csv(file = input.file, header = TRUE, sep = "\t")

str(full.gene.df)
```

    'data.frame':   34254 obs. of  25 variables:
     $ X1                                : chr  "TRINITY_DN37009_c0_g3_i1" "TRINITY_DN37009_c0_g3_i3" "TRINITY_DN37035_c0_g1_i1" "TRINITY_DN37023_c0_g1_i1" ...
     $ X2                                : chr  "sp" "sp" "sp" "sp" ...
     $ X3                                : chr  "Q9BYN0" "Q9BYN0" "P21329" "P07700" ...
     $ X4                                : chr  "SRXN1_HUMAN" "SRXN1_HUMAN" "RTJK_DROFU" "ADRB1_MELGA" ...
     $ X5                                : num  57.6 57.6 28.6 33.3 31.6 ...
     $ X6                                : int  125 125 147 99 133 49 90 179 85 37 ...
     $ X7                                : int  47 47 92 65 87 21 47 105 44 13 ...
     $ X8                                : int  2 2 4 1 3 0 3 7 2 1 ...
     $ X9                                : int  2822 1720 169 47 473 85 15 1 104 102 ...
     $ X10                               : int  2448 1346 579 340 84 231 275 513 343 4 ...
     $ X11                               : int  18 18 613 59 44 374 425 85 64 1097 ...
     $ X12                               : int  136 136 756 157 175 422 509 246 148 1133 ...
     $ X13                               : num  1.11e-39 3.47e-40 4.94e-08 6.59e-09 9.86e-12 ...
     $ X14                               : num  147 147 55.8 55.5 66.6 58.9 49.7 51.6 69.3 46.2 ...
     $ Reviewed                          : chr  "reviewed" "reviewed" NA "reviewed" ...
     $ Entry.Name                        : chr  "SRXN1_HUMAN" "SRXN1_HUMAN" NA "ADRB1_MELGA" ...
     $ Protein.names                     : chr  "Sulfiredoxin-1 (EC 1.8.98.2)" "Sulfiredoxin-1 (EC 1.8.98.2)" NA "Beta-1 adrenergic receptor (Beta-1 adrenoreceptor) (Beta-1 adrenoceptor) (Beta-T)" ...
     $ Gene.Names                        : chr  "SRXN1 C20orf139 SRX SRX1" "SRXN1 C20orf139 SRX SRX1" NA "ADRB1" ...
     $ Organism                          : chr  "Homo sapiens (Human)" "Homo sapiens (Human)" NA "Meleagris gallopavo (Wild turkey)" ...
     $ Length                            : int  137 137 NA 483 483 1056 NA 364 640 1358 ...
     $ Gene.Ontology..biological.process.: chr  "cellular response to oxidative stress [GO:0034599]; response to oxidative stress [GO:0006979]" "cellular response to oxidative stress [GO:0034599]; response to oxidative stress [GO:0006979]" NA "adenylate cyclase-activating adrenergic receptor signaling pathway [GO:0071880]; positive regulation of heart c"| __truncated__ ...
     $ Gene.Ontology..cellular.component.: chr  "cytoplasm [GO:0005737]; cytosol [GO:0005829]; endoplasmic reticulum membrane [GO:0005789]" "cytoplasm [GO:0005737]; cytosol [GO:0005829]; endoplasmic reticulum membrane [GO:0005789]" NA "early endosome [GO:0005769]; membrane [GO:0016020]; plasma membrane [GO:0005886]" ...
     $ Gene.Ontology..GO.                : chr  "cytoplasm [GO:0005737]; cytosol [GO:0005829]; endoplasmic reticulum membrane [GO:0005789]; ATP binding [GO:0005"| __truncated__ "cytoplasm [GO:0005737]; cytosol [GO:0005829]; endoplasmic reticulum membrane [GO:0005789]; ATP binding [GO:0005"| __truncated__ NA "early endosome [GO:0005769]; membrane [GO:0016020]; plasma membrane [GO:0005886]; beta1-adrenergic receptor act"| __truncated__ ...
     $ Gene.Ontology..molecular.function.: chr  "ATP binding [GO:0005524]; oxidoreductase activity, acting on a sulfur group of donors [GO:0016667]; sulfiredoxi"| __truncated__ "ATP binding [GO:0005524]; oxidoreductase activity, acting on a sulfur group of donors [GO:0016667]; sulfiredoxi"| __truncated__ NA "beta1-adrenergic receptor activity [GO:0004940]; identical protein binding [GO:0042802]" ...
     $ Gene.Ontology.IDs                 : chr  "GO:0005524; GO:0005737; GO:0005789; GO:0005829; GO:0006979; GO:0016667; GO:0032542; GO:0034599" "GO:0005524; GO:0005737; GO:0005789; GO:0005829; GO:0006979; GO:0016667; GO:0032542; GO:0034599" NA "GO:0004940; GO:0005769; GO:0005886; GO:0016020; GO:0042802; GO:0045187; GO:0045823; GO:0071880" ...

### 1.1.4 Remove rows with NA, remove whitespace in GO IDs column and keep just gene/GO IDs columns

``` r
# Clean whitespace, filter NA/empty rows, select columns, and split GO terms using column name variables
gene.GO.df <- full.gene.df %>%
  mutate(!!GO.ID.column := str_replace_all(.data[[GO.ID.column]], "\\s*;\\s*", ";")) %>% # Clean up spaces around ";"
  filter(!is.na(.data[[gene.ID.column]]) & !is.na(.data[[GO.ID.column]]) & .data[[GO.ID.column]] != "") %>% 
  select(all_of(c(gene.ID.column, GO.ID.column)))


str(gene.GO.df)
```

    'data.frame':   23192 obs. of  2 variables:
     $ X1               : chr  "TRINITY_DN37009_c0_g3_i1" "TRINITY_DN37009_c0_g3_i3" "TRINITY_DN37023_c0_g1_i1" "TRINITY_DN37023_c1_g1_i1" ...
     $ Gene.Ontology.IDs: chr  "GO:0005524;GO:0005737;GO:0005789;GO:0005829;GO:0006979;GO:0016667;GO:0032542;GO:0034599" "GO:0005524;GO:0005737;GO:0005789;GO:0005829;GO:0006979;GO:0016667;GO:0032542;GO:0034599" "GO:0004940;GO:0005769;GO:0005886;GO:0016020;GO:0042802;GO:0045187;GO:0045823;GO:0071880" "GO:0004940;GO:0005769;GO:0005886;GO:0016020;GO:0042802;GO:0045187;GO:0045823;GO:0071880" ...

### 1.1.5 “Flatten” gene and GO IDs

This flattens the file so all of the GO IDs per gene are separated into
one GO ID per gene per row.

``` r
flat.gene.GO.df <- gene.GO.df %>% separate_rows(!!sym(GO.ID.column), sep = ";")

str(flat.gene.GO.df)
```

    tibble [363,990 × 2] (S3: tbl_df/tbl/data.frame)
     $ X1               : chr [1:363990] "TRINITY_DN37009_c0_g3_i1" "TRINITY_DN37009_c0_g3_i1" "TRINITY_DN37009_c0_g3_i1" "TRINITY_DN37009_c0_g3_i1" ...
     $ Gene.Ontology.IDs: chr [1:363990] "GO:0005524" "GO:0005737" "GO:0005789" "GO:0005829" ...

### 1.1.6 Group

Groups the genes by GO ID (i.e. lists all genes associated with each
unique GO ID)

``` r
grouped.gene.GO.df <- flat.gene.GO.df %>%
  group_by(!!sym(GO.ID.column)) %>%
  summarise(!!gene.ID.column := paste(.data[[gene.ID.column]], collapse = ","))

str(grouped.gene.GO.df)
```

    tibble [15,068 × 2] (S3: tbl_df/tbl/data.frame)
     $ Gene.Ontology.IDs: chr [1:15068] "GO:0000002" "GO:0000009" "GO:0000012" "GO:0000014" ...
     $ X1               : chr [1:15068] "TRINITY_DN12635_c0_g1_i1,TRINITY_DN6578_c1_g1_i2,TRINITY_DN20059_c0_g1_i11,TRINITY_DN20059_c0_g1_i13,TRINITY_DN"| __truncated__ "TRINITY_DN3933_c0_g3_i1" "TRINITY_DN2179_c0_g7_i1,TRINITY_DN2179_c0_g7_i2,TRINITY_DN8923_c0_g2_i1,TRINITY_DN8283_c0_g1_i3,TRINITY_DN8283_"| __truncated__ "TRINITY_DN31365_c0_g1_i1,TRINITY_DN31365_c0_g1_i2,TRINITY_DN53840_c0_g1_i1,TRINITY_DN662_c0_g1_i2,TRINITY_DN149"| __truncated__ ...

### 1.1.7 Map GO IDs to GOslims

The mapping steps were derived from this [bioconductor forum
response](https://support.bioconductor.org/p/128407/#128408)

``` r
# Vector of GO IDs
go_ids <- grouped.gene.GO.df[[GO.ID.column]]

str(go_ids)
```

     chr [1:15068] "GO:0000002" "GO:0000009" "GO:0000012" "GO:0000014" ...

### 1.1.8 Extract GOslims from OBO

Creates new OBO Collection object of just GOslims, based on provided GO
IDs.

``` r
# Create GSEAbase GOCollection using `go_ids`
myCollection <- GOCollection(go_ids)

# Retrieve GOslims from GO OBO file set
slim <- getOBOCollection(gseabase_files)

str(slim)
```

    Formal class 'OBOCollection' [package "GSEABase"] with 7 slots
      ..@ .stanza     :'data.frame':    153 obs. of  1 variable:
      .. ..$ value: chr [1:153] "Root" "Term" "Term" "Term" ...
      ..@ .subset     :'data.frame':    22 obs. of  1 variable:
      .. ..$ value: chr [1:22] "Rhea list of ChEBI terms representing the major species at pH 7.3." "Term not to be used for direct annotation" "Terms planned for obsoletion" "AGR slim" ...
      ..@ .kv         :'data.frame':    2075 obs. of  3 variables:
      .. ..$ stanza_id: chr [1:2075] ".__Root__" ".__Root__" ".__Root__" ".__Root__" ...
      .. ..$ key      : chr [1:2075] "format-version" "data-version" "synonymtypedef" "synonymtypedef" ...
      .. ..$ value    : chr [1:2075] "1.2" "go/2024-09-08/subsets/goslim_generic.owl" "syngo_official_label \"label approved by the SynGO project\"" "systematic_synonym \"Systematic synonym\" EXACT" ...
      ..@ evidenceCode: chr [1:26] "EXP" "IDA" "IPI" "IMP" ...
      ..@ ontology    : chr NA
      ..@ ids         : chr [1:141] "GO:0000228" "GO:0000278" "GO:0000910" "GO:0001618" ...
      ..@ type        : chr "OBO"

### 1.1.9 Get Biological Process (BP) GOslims associated with provided GO IDs.

``` r
# Retrieve Biological Process (BP) GOslims
slimdf <- goSlim(myCollection, slim, "BP", verbose)
str(slimdf)
```

    'data.frame':   72 obs. of  3 variables:
     $ Count  : int  93 22 18 554 81 107 31 5 103 58 ...
     $ Percent: num  0.925 0.219 0.179 5.508 0.805 ...
     $ Term   : chr  "mitotic cell cycle" "cytokinesis" "cytoplasmic translation" "immune system process" ...

### 1.1.10 Map GO to GOslims

Performs mapping of of GOIDs to GOslims

Returns:

- GOslim IDs (as rownames)
- GOslim terms
- Counts of GO IDs matching to corresponding GOslim
- Percentage of GO IDs matching to corresponding GOslim
- GOIDs mapped to corresponding GOslim, in a semi-colon delimited format

``` r
# List of GOslims and all GO IDs from `go_ids`
gomap <- as.list(GOBPOFFSPRING[rownames(slimdf)])

# Maps `go_ids` to matching GOslims
mapped <- lapply(gomap, intersect, ids(myCollection))

# Append all mapped GO IDs to `slimdf`
# `sapply` needed to apply paste() to create semi-colon delimited values
slimdf$GO.IDs <- sapply(lapply(gomap, intersect, ids(myCollection)), paste, collapse=";")

# Remove "character(0) string from "GO.IDs" column
slimdf$GO.IDs[slimdf$GO.IDs == "character(0)"] <- ""

# Add self-matching GOIDs to "GO.IDs" column, if not present
for (go_id in go_ids) {
  # Check if the go_id is present in the row names
  if (go_id %in% rownames(slimdf)) {
    # Check if the go_id is not present in the GO.IDs column
    # Also removes white space "trimws()" and converts all to upper case to handle
    # any weird, "invisible" formatting issues.
    if (!go_id %in% trimws(toupper(strsplit(slimdf[go_id, "GO.IDs"], ";")[[1]]))) {
      # Append the go_id to the GO.IDs column with a semi-colon separator
      if (length(slimdf$GO.IDs) > 0 && nchar(slimdf$GO.IDs[nrow(slimdf)]) > 0) {
        slimdf[go_id, "GO.IDs"] <- paste0(slimdf[go_id, "GO.IDs"], "; ", go_id)
      } else {
        slimdf[go_id, "GO.IDs"] <- go_id
      }
    }
  }
}

str(slimdf)
```

    'data.frame':   72 obs. of  4 variables:
     $ Count  : int  93 22 18 554 81 107 31 5 103 58 ...
     $ Percent: num  0.925 0.219 0.179 5.508 0.805 ...
     $ Term   : chr  "mitotic cell cycle" "cytokinesis" "cytoplasmic translation" "immune system process" ...
     $ GO.IDs : chr  "GO:0000022;GO:0000070;GO:0000082;GO:0000086;GO:0000132;GO:0000281;GO:0006977;GO:0007052;GO:0007076;GO:0007079;G"| __truncated__ "GO:0000281;GO:0000911;GO:0000915;GO:0007112;GO:0031991;GO:0032465;GO:0032466;GO:0032467;GO:0032506;GO:0033206;G"| __truncated__ "GO:0001731;GO:0001732;GO:0002183;GO:0002184;GO:0002188;GO:0002190;GO:0002191;GO:0002192;GO:0140018;GO:0140708;G"| __truncated__ "GO:0001768;GO:0001771;GO:0001773;GO:0001774;GO:0001776;GO:0001779;GO:0001780;GO:0001782;GO:0001798;GO:0001805;G"| __truncated__ ...

### 1.1.11 Flatten GOslims

``` r
# "Flatten" file so each row is single GO ID with corresponding GOslim
# rownames_to_column needed to retain row name info
slimdf_separated <- as.data.frame(slimdf %>%
  rownames_to_column('GOslim') %>%
  separate_rows(GO.IDs, sep = ";"))

# Group by unique GO ID
grouped_slimdf <- slimdf_separated %>%
  filter(!is.na(GO.IDs) & GO.IDs != "") %>%
  group_by(GO.IDs) %>%
  summarize(GOslim = paste(GOslim, collapse = ";"),
            Term = paste(Term, collapse = ";"))


str(grouped_slimdf)
```

    tibble [6,942 × 3] (S3: tbl_df/tbl/data.frame)
     $ GO.IDs: chr [1:6942] " GO:0000278" " GO:0002181" " GO:0002376" " GO:0003014" ...
     $ GOslim: chr [1:6942] "GO:0000278" "GO:0002181" "GO:0002376" "GO:0003014" ...
     $ Term  : chr [1:6942] "mitotic cell cycle" "cytoplasmic translation" "immune system process" "renal system process" ...

## 1.2 Write slimdf to file

### 1.2.1 Sort and select slmidf

Sorts GOslims by `Count`, in descending order and then selects just the
`Term` and `Count` columns.

``` r
slimdf.sorted <- slimdf %>% arrange(desc(Count))

slim.count.df <- slimdf.sorted %>% 
  select(Term, Count)

str(slim.count.df)
```

    'data.frame':   72 obs. of  2 variables:
     $ Term : chr  "anatomical structure development" "cell differentiation" "signaling" "immune system process" ...
     $ Count: int  2091 932 908 554 352 306 278 275 264 246 ...

### 1.2.2 Write formatted slim.count.df to file

``` r
write_tsv(slim.count.df, file = "../output/01-GO-to-GOslim/counts.GOID-per-GOslim_term.tab")
```

## 1.3 Count unique Biological Process GO IDs

``` r
# Flatten the list and count total GO IDs
total_go_ids <- length(unlist(gomap))

# Display the total count
total_go_ids
```

    [1] 26534

``` r
# Unlist to extract all GO IDs, then find unique ones and count them
unique_go_ids <- unique(unlist(gomap))
total_unique_ids <- length(unique_go_ids)

# Display the total count of unique GO IDs
total_unique_ids
```

    [1] 17913

Total starting BP GO IDs: 26534

Total unique BP GO IDs: 17913

# 2 REFERENCES

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-thegene2023" class="csl-entry">

Aleksander, Suzi A, James Balhoff, Seth Carbon, J Michael Cherry, Harold
J Drabkin, Dustin Ebert, Marc Feuermann, et al. 2023. “The Gene Ontology
Knowledgebase in 2023.” Edited by A Baryshnikova. *GENETICS* 224 (1).
<https://doi.org/10.1093/genetics/iyad031>.

</div>

<div id="ref-ashburner2000" class="csl-entry">

Ashburner, Michael, Catherine A. Ball, Judith A. Blake, David Botstein,
Heather Butler, J. Michael Cherry, Allan P. Davis, et al. 2000. “Gene
Ontology: Tool for the Unification of Biology.” *Nature Genetics* 25
(1): 25–29. <https://doi.org/10.1038/75556>.

</div>

<div id="ref-carlson2017" class="csl-entry">

Carlson, Marc. 2017. “GO.db.” <https://doi.org/10.18129/B9.BIOC.GO.DB>.

</div>

<div id="ref-martinmorgan2017" class="csl-entry">

Martin Morgan, Seth Falcon. 2017. “GSEABase.”
<https://doi.org/10.18129/B9.BIOC.GSEABASE>.

</div>

<div id="ref-wickham2019" class="csl-entry">

Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy
McGowan, Romain François, Garrett Grolemund, et al. 2019. “Welcome to
the Tidyverse.” *Journal of Open Source Software* 4 (43): 1686.
<https://doi.org/10.21105/joss.01686>.

</div>

</div>
