`20230328-pgen-gene_annotation-update/scripts`

---

- [01-uniprot_parsing-bash_variables.txt](https://github.com/RobertsLab/code/blob/master/r_projects/sam/20230328-pgen-gene_annotation-update/scripts/01-uniprot_parsing-bash_variables.txt): Text file containing `bash` variables to use throughout the [01-uniprot_parsing.Rmd](https://github.com/RobertsLab/code/blob/master/r_projects/sam/20230328-pgen-gene_annotation-update/scripts/01-uniprot_parsing.Rmd) R Markdown script.

- [01-uniprot_parsing.Rmd](https://github.com/RobertsLab/code/blob/master/r_projects/sam/20230328-pgen-gene_annotation-update/scripts/01-uniprot_parsing.Rmd): R Markdown file which parses UniProt batch retrieval file, generated on [20220419 by Sam White](https://robertslab.github.io/sams-notebook/2022/04/19/Data-Wrangling-Create-Primary-P.generosa-Genome-Annotation-File.html) (notebook entry), for UniProt accessions, gene IDs, gene names, and gene ontology (GO) IDs.

- [02-goslim-mapping.Rmd](https://github.com/RobertsLab/code/blob/master/r_projects/sam/20230328-pgen-gene_annotation-update/scripts/02-goslim-mapping.Rmd): R Markdown file which maps gene ontology (GO) IDs identified in [01-uniprot_parsing.Rmd](https://github.com/RobertsLab/code/blob/master/r_projects/sam/20230328-pgen-gene_annotation-update/scripts/01-uniprot_parsing.Rmd) to Biological Process (BP) GOslims. Then, groups/joins results by [_Panopea generosa_ (Pacific geoduck)](http://en.wikipedia.org/wiki/Geoduck) genes.