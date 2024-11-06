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

- [pyc_go.tsv](https://gannet.fish.washington.edu/seashell/bu-github/project-pycno-multispecies-2023/output/06-Go-3-species/pyc_go.tsv)

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
input.file <- "https://gannet.fish.washington.edu/seashell/bu-github/project-pycno-multispecies-2023/output/06-Go-3-species/pyc_go.tsv"


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

    'data.frame':   29145 obs. of  25 variables:
     $ X1                                : chr  "TRINITY_DN37046_c0_g1_i2" "TRINITY_DN37065_c0_g1_i1" "TRINITY_DN37015_c0_g1_i1" "TRINITY_DN45769_c0_g1_i1" ...
     $ X2                                : chr  "sp" "sp" "sp" "sp" ...
     $ X3                                : chr  "Q32NR4" "Q95NT6" "Q9Z0Z4" "Q09575" ...
     $ X4                                : chr  "TTC29_XENLA" "MTH2_DROYA" "HEPH_MOUSE" "YRD6_CAEEL" ...
     $ X5                                : num  36.4 26 47.5 37.2 51.7 ...
     $ X6                                : int  439 231 80 121 58 60 35 88 78 88 ...
     $ X7                                : int  269 155 42 75 28 35 17 59 42 45 ...
     $ X8                                : int  5 8 0 1 0 1 0 1 1 3 ...
     $ X9                                : int  502 688 2 361 12 217 55 86 274 266 ...
     $ X10                               : int  1818 29 241 2 185 38 159 349 41 12 ...
     $ X11                               : int  19 224 38 776 1146 445 297 456 437 34 ...
     $ X12                               : int  447 449 117 896 1203 503 331 542 512 113 ...
     $ X13                               : num  6.37e-79 2.22e-13 3.37e-21 4.39e-25 1.02e-14 ...
     $ X14                               : num  263 75.5 89 102 70.1 55.8 45.8 50.1 62.8 47.8 ...
     $ Reviewed                          : chr  NA NA "reviewed" NA ...
     $ Entry.Name                        : chr  NA NA "HEPH_MOUSE" NA ...
     $ Protein.names                     : chr  NA NA "Hephaestin (Hp) (EC 1.16.3.1)" NA ...
     $ Gene.Names                        : chr  NA NA "Heph Kiaa0698" NA ...
     $ Organism                          : chr  NA NA "Mus musculus (Mouse)" NA ...
     $ Length                            : int  NA NA 1157 NA 2812 NA NA NA 583 NA ...
     $ Gene.Ontology..biological.process.: chr  NA NA "intestinal iron absorption [GO:0160179]; iron ion transport [GO:0006826]; multicellular organismal-level iron i"| __truncated__ NA ...
     $ Gene.Ontology..cellular.component.: chr  NA NA "basolateral plasma membrane [GO:0016323]; extracellular region [GO:0005576]; perinuclear region of cytoplasm [G"| __truncated__ NA ...
     $ Gene.Ontology..GO.                : chr  NA NA "basolateral plasma membrane [GO:0016323]; extracellular region [GO:0005576]; perinuclear region of cytoplasm [G"| __truncated__ NA ...
     $ Gene.Ontology..molecular.function.: chr  NA NA "copper ion binding [GO:0005507]; ferrous iron binding [GO:0008198]; ferroxidase activity [GO:0004322]; oxidored"| __truncated__ NA ...
     $ Gene.Ontology.IDs                 : chr  NA NA "GO:0004322; GO:0005507; GO:0005576; GO:0005886; GO:0006826; GO:0008198; GO:0016323; GO:0016491; GO:0048471; GO:"| __truncated__ NA ...

### 1.1.4 Remove rows with NA, remove whitespace in GO IDs column and keep just gene/GO IDs columns

``` r
# Clean whitespace, filter NA/empty rows, select columns, and split GO terms using column name variables
gene.GO.df <- full.gene.df %>%
  mutate(!!GO.ID.column := str_replace_all(.data[[GO.ID.column]], "\\s*;\\s*", ";")) %>% # Clean up spaces around ";"
  filter(!is.na(.data[[gene.ID.column]]) & !is.na(.data[[GO.ID.column]]) & .data[[GO.ID.column]] != "") %>% 
  select(all_of(c(gene.ID.column, GO.ID.column)))


str(gene.GO.df)
```

    'data.frame':   19842 obs. of  2 variables:
     $ X1               : chr  "TRINITY_DN37015_c0_g1_i1" "TRINITY_DN45794_c0_g1_i1" "TRINITY_DN45764_c0_g1_i1" "TRINITY_DN48070_c0_g1_i1" ...
     $ Gene.Ontology.IDs: chr  "GO:0004322;GO:0005507;GO:0005576;GO:0005886;GO:0006826;GO:0008198;GO:0016323;GO:0016491;GO:0048471;GO:0060586;G"| __truncated__ "GO:0005615;GO:0005886;GO:0007339;GO:0031012;GO:0098609" "GO:0000139;GO:0001736;GO:0001737;GO:0004674;GO:0005615;GO:0005737;GO:0005886;GO:0006468;GO:0007219;GO:0007267;G"| __truncated__ "GO:0005102;GO:0005739;GO:0005777;GO:0005778;GO:0006633;GO:0019166;GO:0033306;GO:0043231" ...

### 1.1.5 “Flatten” gene and GO IDs

This flattens the file so all of the GO IDs per gene are separated into
one GO ID per gene per row.

``` r
flat.gene.GO.df <- gene.GO.df %>% separate_rows(!!sym(GO.ID.column), sep = ";")

str(flat.gene.GO.df)
```

    tibble [319,980 × 2] (S3: tbl_df/tbl/data.frame)
     $ X1               : chr [1:319980] "TRINITY_DN37015_c0_g1_i1" "TRINITY_DN37015_c0_g1_i1" "TRINITY_DN37015_c0_g1_i1" "TRINITY_DN37015_c0_g1_i1" ...
     $ Gene.Ontology.IDs: chr [1:319980] "GO:0004322" "GO:0005507" "GO:0005576" "GO:0005886" ...

### 1.1.6 Group

Groups the genes by GO ID (i.e. lists all genes associated with each
unique GO ID)

``` r
grouped.gene.GO.df <- flat.gene.GO.df %>%
  group_by(!!sym(GO.ID.column)) %>%
  summarise(!!gene.ID.column := paste(.data[[gene.ID.column]], collapse = ","))

str(grouped.gene.GO.df)
```

    tibble [14,922 × 2] (S3: tbl_df/tbl/data.frame)
     $ Gene.Ontology.IDs: chr [1:14922] "GO:0000002" "GO:0000009" "GO:0000011" "GO:0000012" ...
     $ X1               : chr [1:14922] "TRINITY_DN624_c0_g1_i1,TRINITY_DN624_c0_g1_i2,TRINITY_DN678_c0_g1_i1,TRINITY_DN678_c0_g1_i2,TRINITY_DN678_c0_g1"| __truncated__ "TRINITY_DN10071_c0_g4_i1" "TRINITY_DN628_c0_g2_i1,TRINITY_DN9588_c0_g1_i1,TRINITY_DN9588_c0_g1_i5" "TRINITY_DN3105_c1_g1_i3,TRINITY_DN4630_c0_g1_i2,TRINITY_DN4630_c0_g1_i4,TRINITY_DN4630_c0_g1_i5,TRINITY_DN8346_"| __truncated__ ...

### 1.1.7 Map GO IDs to GOslims

The mapping steps were derived from this [bioconductor forum
response](https://support.bioconductor.org/p/128407/#128408)

``` r
# Vector of GO IDs
go_ids <- grouped.gene.GO.df[[GO.ID.column]]

str(go_ids)
```

     chr [1:14922] "GO:0000002" "GO:0000009" "GO:0000011" "GO:0000012" ...

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
     $ Count  : int  102 23 17 533 83 107 28 5 96 60 ...
     $ Percent: num  1.016 0.229 0.169 5.308 0.827 ...
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
     $ Count  : int  102 23 17 533 83 107 28 5 96 60 ...
     $ Percent: num  1.016 0.229 0.169 5.308 0.827 ...
     $ Term   : chr  "mitotic cell cycle" "cytokinesis" "cytoplasmic translation" "immune system process" ...
     $ GO.IDs : chr  "GO:0000022;GO:0000070;GO:0000082;GO:0000086;GO:0000132;GO:0000281;GO:0000282;GO:0001100;GO:0006977;GO:0007052;G"| __truncated__ "GO:0000281;GO:0000282;GO:0000911;GO:0000915;GO:0000917;GO:0007111;GO:0032465;GO:0032466;GO:0032467;GO:0032506;G"| __truncated__ "GO:0001731;GO:0001732;GO:0002182;GO:0002183;GO:0002190;GO:0002191;GO:0002192;GO:0140018;GO:0140708;GO:1903679;G"| __truncated__ "GO:0001771;GO:0001773;GO:0001774;GO:0001776;GO:0001777;GO:0001779;GO:0001780;GO:0001782;GO:0001805;GO:0001812;G"| __truncated__ ...

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

    tibble [6,903 × 3] (S3: tbl_df/tbl/data.frame)
     $ GO.IDs: chr [1:6903] " GO:0000278" " GO:0002181" " GO:0002376" " GO:0003014" ...
     $ GOslim: chr [1:6903] "GO:0000278" "GO:0002181" "GO:0002376" "GO:0003014" ...
     $ Term  : chr [1:6903] "mitotic cell cycle" "cytoplasmic translation" "immune system process" "renal system process" ...

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
     $ Count: int  2120 941 925 533 331 316 268 267 253 244 ...

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
