20240916-snam-c1q-primer-design
================
Sam White
2024-09-16

- <a href="#1-create-bash-variables-file"
  id="toc-1-create-bash-variables-file">1 CREATE BASH VARIABLES FILE</a>
- <a href="#2-download-ncbi-genome-files"
  id="toc-2-download-ncbi-genome-files">2 DOWNLOAD NCBI GENOME FILES</a>
  - <a href="#21-download-the-actual-files"
    id="toc-21-download-the-actual-files">2.1 Download the actual files</a>
  - <a href="#22-check-md5-checkums" id="toc-22-check-md5-checkums">2.2
    Check MD5 Checkums</a>
  - <a href="#23-decompress-ncbi-files"
    id="toc-23-decompress-ncbi-files">2.3 Decompress NCBI files</a>
- <a href="#3-extract-c1q-gene-sequence"
  id="toc-3-extract-c1q-gene-sequence">3 EXTRACT C1Q GENE SEQUENCE</a>
  - <a href="#31-peek-at-gff" id="toc-31-peek-at-gff">3.1 Peek at GFF</a>
  - <a href="#32-format-region-and-extract-c1q-sequence-as-fasta"
    id="toc-32-format-region-and-extract-c1q-sequence-as-fasta">3.2 Format
    region and extract C1q sequence as FastA.</a>
- <a href="#4-extract-c1q-with-53-buffer-regions"
  id="toc-4-extract-c1q-with-53-buffer-regions">4 EXTRACT C1Q WITH 5’/3’
  BUFFER REGIONS</a>
- <a href="#5-primer-design-using-primer3"
  id="toc-5-primer-design-using-primer3">5 PRIMER DESIGN USING PRIMER3</a>
  - <a href="#51-design-primers" id="toc-51-design-primers">5.1 Design
    primers</a>
  - <a href="#52-amplicon-length" id="toc-52-amplicon-length">5.2 Amplicon
    length</a>
- <a href="#6-split-genome" id="toc-6-split-genome">6 SPLIT GENOME</a>
  - <a
    href="#61-split-mulit-fasta-file-in-to-individual-fasta-files-with-pyfaidx"
    id="toc-61-split-mulit-fasta-file-in-to-individual-fasta-files-with-pyfaidx">6.1
    Split mulit-FastA file in to individual FastA files with PyFaidx</a>
- <a href="#7-primer-search-with-emboss-primersearch"
  id="toc-7-primer-search-with-emboss-primersearch">7 PRIMER SEARCH WITH
  EMBOSS PRIMERSEARCH</a>
  - <a href="#71-create-emboss-primersearch-primers-file"
    id="toc-71-create-emboss-primersearch-primers-file">7.1 Create EMBOSS
    PrimerSearch Primers File</a>
  - <a href="#72-run-emboss-primersearch"
    id="toc-72-run-emboss-primersearch">7.2 Run EMBOSS PrimerSearch</a>
  - <a href="#73-check-primer-matches" id="toc-73-check-primer-matches">7.3
    Check primer matches</a>
- <a href="#8-summary" id="toc-8-summary">8 SUMMARY</a>
- <a href="#9-citations" id="toc-9-citations">9 CITATIONS</a>

This notebook describes using
[Primer3](https://github.com/primer3-org/primer3) ([Untergasser et al.
2012](#ref-untergasser2012); [Koressaar and Remm
2007](#ref-koressaar2007)) to reproducibly design primers for sequencing
of the [lake trout
(*S.namaycush*)](https://en.wikipedia.org/wiki/Lake_trout) (Wikipedia)
C1q gene
([`LOC120027825`](https://github.com/RobertsLab/resources/issues/1941#issuecomment-2246548960)
(GitHub Issue)). Additionally, this is primarily to confirm that we can
successfully amplify the C1q gene, since the [bisulfite PCR was
unsuccessful](https://robertslab.github.io/sams-notebook/posts/2024/2024-09-11-Bisulfite-PCR---Second-Primer-Annealing-Gradient-Test-with-Lake-Trout-Bisulfite-treated-DNA/)
(Notebook).

This process also utilizes [pyfaidx](https://github.com/mdshw5/pyfaidx)
([Shirley et al. 2015](#ref-shirley2015)).

[EMBOSS
PrimerSearch](https://emboss.sourceforge.net/apps/cvs/emboss/apps/primersearch.html)
([Rice, Longden, and Bleasby 2000](#ref-rice2000)) will be utilized to
assess primer specificity across the genome.

The notebook entry was rendered directly from the R Markdown file,
[`20240916-snam-c1q-primer-design.Rmd`](https://github.com/RobertsLab/code/blob/ea893a0e9607662b85278bb620f042994efbf374/r_projects/sam/20240916-snam-c1q-primer-design/code/20240916-snam-c1q-primer-design.Rmd)
(GitHub - commit ea893a0e9607662b85278bb620f042994efbf374).

# 1 CREATE BASH VARIABLES FILE

This allows usage of Bash variables across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# DATA DIRECTORIES"
echo 'export data_dir="../data/S_namaycush/genomes"'
echo 'export output_top="../output"'
echo 'export genome_fasta_splits_dir="${data_dir}/fasta_splits"'
echo ""

echo "# SEQUENCE"
echo 'export sequence_ID="LOC120027825"'
echo ""

echo "# SEQUENCE REGIONS"
echo 'export left_buffer="500"'
echo 'export right_buffer="500"'
echo ""

echo "# INPUT FILES"
echo 'export genome_fasta="GCF_016432855.1_SaNama_1.0_genomic.fna"'
echo 'export genome_gff="GCF_016432855.1_SaNama_1.0_genomic.gff"'
echo 'export ncbi_gff_gz="GCF_016432855.1_SaNama_1.0_genomic.gff.gz"'
echo 'export ncbi_fasta_gz="GCF_016432855.1_SaNama_1.0_genomic.fna.gz"'
echo 'export ncbi_md5sums="md5checksums.txt"'
echo 'export ncbi_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/432/855/GCF_016432855.1_SaNama_1.0/"'
echo ""

echo "# OUTPUT FILES"
echo 'export c1q_fasta="${sequence_ID}.fasta"'
echo 'export c1q_buffer_fasta="${sequence_ID}-${left_buffer}bp.left-${right_buffer}bp.right.fasta"'
echo 'export c1q_faidx_region_file="${sequence_ID}"-region.txt'
echo 'export c1q_buffer_faidx_region_file="${sequence_ID}-${left_buffer}bp.left-${right_buffer}bp.right-region.txt"'
echo ""

echo "# SET CPUS"
echo 'export threads=40'
echo ""

echo "# PROGRAMS"
echo 'export pyfaidx=/home/shared/pyfaidx-0.8.1.1'
echo 'export primer3_dir="/home/shared/primer3-2.6.1/src"'
echo 'export primer3="${primer3_dir}/primer3_core"'
echo 'export primer3_config="${primer3_dir}/primer3_config"'
echo 'export primersearch="/home/shared/EMBOSS-6.6.0/emboss/primersearch"'



} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # DATA DIRECTORIES
    export data_dir="../data/S_namaycush/genomes"
    export output_top="../output"
    export genome_fasta_splits_dir="${data_dir}/fasta_splits"

    # SEQUENCE
    export sequence_ID="LOC120027825"

    # SEQUENCE REGIONS
    export left_buffer="500"
    export right_buffer="500"

    # INPUT FILES
    export genome_fasta="GCF_016432855.1_SaNama_1.0_genomic.fna"
    export genome_gff="GCF_016432855.1_SaNama_1.0_genomic.gff"
    export ncbi_gff_gz="GCF_016432855.1_SaNama_1.0_genomic.gff.gz"
    export ncbi_fasta_gz="GCF_016432855.1_SaNama_1.0_genomic.fna.gz"
    export ncbi_md5sums="md5checksums.txt"
    export ncbi_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/432/855/GCF_016432855.1_SaNama_1.0/"

    # OUTPUT FILES
    export c1q_fasta="${sequence_ID}.fasta"
    export c1q_buffer_fasta="${sequence_ID}-${left_buffer}bp.left-${right_buffer}bp.right.fasta"
    export c1q_faidx_region_file="${sequence_ID}"-region.txt
    export c1q_buffer_faidx_region_file="${sequence_ID}-${left_buffer}bp.left-${right_buffer}bp.right-region.txt"

    # SET CPUS
    export threads=40

    # PROGRAMS
    export pyfaidx=/home/shared/pyfaidx-0.8.1.1
    export primer3_dir="/home/shared/primer3-2.6.1/src"
    export primer3="${primer3_dir}/primer3_core"
    export primer3_config="${primer3_dir}/primer3_config"
    export primersearch="/home/shared/EMBOSS-6.6.0/emboss/primersearch"

# 2 DOWNLOAD NCBI GENOME FILES

## 2.1 Download the actual files

``` bash
# Load bash variables into memory
source .bashvars

for file in ${ncbi_gff_gz} ${ncbi_fasta_gz} ${ncbi_md5sums}
do
  wget \
  --no-check-certificate \
  --continue \
  --quiet \
  --directory-prefix=${data_dir} \
  ${ncbi_url}${file}
done

ls -lh "${data_dir}"
```

    total 3.3G
    drwxr-xr-x 2 sam sam 164K Sep 16 16:44 fasta_splits
    -rw-r--r-- 1 sam sam 2.3G Jan 13  2021 GCF_016432855.1_SaNama_1.0_genomic.fna
    -rw-r--r-- 1 sam sam 154K Sep 16 13:22 GCF_016432855.1_SaNama_1.0_genomic.fna.fai
    -rw-r--r-- 1 sam sam 642M Jan 13  2021 GCF_016432855.1_SaNama_1.0_genomic.fna.gz
    -rw-r--r-- 1 sam sam 373M Sep 28  2022 GCF_016432855.1_SaNama_1.0_genomic.gff
    -rw-r--r-- 1 sam sam  18M Sep 28  2022 GCF_016432855.1_SaNama_1.0_genomic.gff.gz
    -rw-r--r-- 1 sam sam  33K Sep 16 01:22 md5checksums.txt

## 2.2 Check MD5 Checkums

``` bash
# Load bash variables into memory
source .bashvars

cd "${data_dir}"

for file in *.gz
do
  grep "${file}" ${ncbi_md5sums} | md5sum -c -
done
```

    ./GCF_016432855.1_SaNama_1.0_genomic.fna.gz: OK
    ./GCF_016432855.1_SaNama_1.0_genomic.gff.gz: OK

## 2.3 Decompress NCBI files

``` bash
# Load bash variables into memory
source .bashvars

cd "${data_dir}"

for file in *.gz
do
  gunzip "${file}"
done

ls -lh
```

    gzip: GCF_016432855.1_SaNama_1.0_genomic.fna already exists;    not overwritten
    gzip: GCF_016432855.1_SaNama_1.0_genomic.gff already exists;    not overwritten
    total 3.3G
    drwxr-xr-x 2 sam sam 164K Sep 16 16:44 fasta_splits
    -rw-r--r-- 1 sam sam 2.3G Jan 13  2021 GCF_016432855.1_SaNama_1.0_genomic.fna
    -rw-r--r-- 1 sam sam 154K Sep 16 13:22 GCF_016432855.1_SaNama_1.0_genomic.fna.fai
    -rw-r--r-- 1 sam sam 642M Jan 13  2021 GCF_016432855.1_SaNama_1.0_genomic.fna.gz
    -rw-r--r-- 1 sam sam 373M Sep 28  2022 GCF_016432855.1_SaNama_1.0_genomic.gff
    -rw-r--r-- 1 sam sam  18M Sep 28  2022 GCF_016432855.1_SaNama_1.0_genomic.gff.gz
    -rw-r--r-- 1 sam sam  33K Sep 16 01:22 md5checksums.txt

# 3 EXTRACT C1Q GENE SEQUENCE

## 3.1 Peek at GFF

``` bash
# Load bash variables into memory
source .bashvars

awk '$3=="gene"' "${data_dir}/${genome_gff}" | grep "${sequence_ID}" | column -t
```

    NC_052339.1  Gnomon  gene  19244556  19248013  .  -  .  ID=gene-LOC120027825;Dbxref=GeneID:120027825;Name=LOC120027825;gbkey=Gene;gene=LOC120027825;gene_biotype=protein_coding

## 3.2 Format region and extract C1q sequence as FastA.

``` bash
# Load bash variables into memory
source .bashvars

# Format region for use with pyfaidx
region=$(awk '$3=="gene"' "${data_dir}/${genome_gff}" | grep "${sequence_ID}" | awk '{print $1":"$4"-"$5}')

echo "${region}"
echo ""

# Extract region
${pyfaidx} "${data_dir}/${genome_fasta}" "${region}" | tee "${output_top}/${c1q_fasta}"
```

    NC_052339.1:19244556-19248013

    >NC_052339.1:19244556-19248013
    TGGGATGTGTAATGGTCATATTTATTACCATACCTCTCCACTTCTTCCTCTTTCCctttacaatttttattttttttaaa
    taataaagtGATGATTTATTTATTGGTTAAAGAGGAAAATCCTCACATCACAGAATGGTGCTGTCCTCACACTGGGAAGA
    GCAGGAAGCCACTGAAGGTGTTGTGGTTATCATGGCTATCATGGAGACCCTGGTCTTCAGGGAGACGAAGGTAGACCACA
    TCCTCCTTCTCCAGCTCTAGTGTCAACGCATTAGATATGTACTGCCAACCCCCAACTGTATTTCTTTCTACATTATATAG
    AACTCTTTGATTATTGTGAAACATCATGATGCCCATGAGTTGTGGGTGGCGGGATGACATGGCCGTGAATCTGAAGAAGT
    AGACTCCTTTCACTGATGCTGTGAAGATGCCTGCATCAGAGGAGAAATAGAGTGATATCAGCACTGTCTGTCTCATGTTG
    CTGTAGCTGGAATGTTTTtaacaaataaactcagcaaaaaaagaaatgtcccttttcaggaccctgtctttcaaagataa
    ttcgtaaaaatccaaataacttcacagatcttcattttaaagggtttaaacactgtttcccatgcttcttcaatgaacca
    taaacaattaatgaacatgcacctgtggaacggtcgttaagacactaacagcttacagacggtaggcagttAGGCcaaag
    ttatgaaaacttaggacactaaagaggcctttctactgacactgaaaaacaccaaaagaaagatgtccagggtccctgct
    catctgcaggaacgtgccttagacatgctgcaaggaggcatgaggactacagatgtggccagggcaataaattgcaatgt
    ccgtactgtgagacacctaacacagcgctacagggagacaggacggacagctgatcgtcctcgcagtgacagaccacgtg
    taacaacacctgcacaggatcggtacatccgaacatcacacctgcgggacaggtaccggatggcaacaacaactgcccga
    gttacaccaggaacgcacaatccctccatcagtgctcagactgtctgcaataggctaagagaggttggactgagggcttg
    taggcctgttgtcaGGCAGGTCATCACCGGCAActacgtcgcctatgggcacaaacccaccgtcactggaccagacagcc
    gcggttttgtctcaccaggggtgatggtcggaattgcgtttatcgttgaaggaatgagcgttacactgaggcctgtactc
    tggagcgggatcgatttggaggtggagggtccgtcatggtctggggcggtgtatcacagcatcatcggactgagcttgtt
    gtctttgcaggaaatctcaacgctgtgcgttacagggaagacatcctcctccctcatatggtacccttcctgcaggctca
    tcctgacatgaccctccagcatgacaatgccaccagccatactgctcatgttgtgcgtgatttcctgcaagacaggaatg
    tcagtgttctgccatggccagcgaagagcccagatctcaatcccattgagcacgtctgggaccttttggatcggagggtg
    agagctagggccattcgccacagaaatgtctgggaacttgcaggtgccttggtggaagagtggggtaacatctcacagca
    agaactggcaaatctggtgcagtccatgaggaggagatgcactgcagtacttaatgcagctggtggccacaccagatact
    gactgttacttttgattttgaccccacctttgttcagggacacattattccatttctgttagtcacatgtctgtggaact
    tgttcagtttatgtctcagttgttgaatattgttatgttcatacaaatatttacacatgttaagttgactgaaaataaac
    gcagttgacagtgagaggacgtttcttttaaCAGTGAGAGGAGTTTAGTATAATTAAGTCAGGGACTAGGTAACACAAAT
    TACCACCAAACACAATAAGTATCAAACTAAACATTGTTATAACTAAAATAAAATCAGGGTGACAGATTAAGTAACATAAA
    ATATGTCATTAAAAAAAATGACCGTAgacatacagtactagtcaaaagtttggacatacttagaggtgtgggtttttctt
    tatttttactattttctacattgtagcataatagtgaagacataaactatgaaataacaaatatggaatcatgtagtaac
    caaaatagtatttaaaaaatctaaatatatttcagatttttcaaagtagccaccctttgcctcgatgacagctttgcact
    ctcttggcattctctcttcacctggaatgcttttccaagagtgtgcaaaactgtcatcaaggcaaagggtgactacttta
    aagaatgtcaaatataaaatatattttgatttgttttatactttttgggttactacatgattccatatgtgttatttaat
    agttttgaagtcttcactattattctacaatgtagaaaatagtaaaataaagaatccctggaataagtaggtgtgtccaa
    acttttgactggtactgtacgtgttTGGTATCATATAGCTAGTAGTTACAAAAGTAAACTCATTTAAAAGCTTGAAACGT
    AGAATTTGTCCCTCAGAATTATGTCCACTCTGAAAATCTAGTGTCTTATCATGATTTTCAAAGTAAAATATAAAACAAGT
    TAGTGAAGTAAAATTATATACAGAGATTACTGAGAAACATGCAGCAAAATCATCCCAGAAGAATTACATTAAATTGCATT
    CATAAAAAATTATGATGAGATGACTGACTCAATATTTGGCATTGTTGGATGAATGTTATGCAACATTCAAATGGATTTCC
    CTCTTGTTAGTTAGTACCTGTAATTGGGTTGTAGGCCTTGCCGATGTTGGTGATGACTCTACTGAAGATCAGTGTAGTAG
    CAGTACTGAATGGACCTACGTTTCCAGAGTCAGTCAAACCAGCAGAGAAGGCCACAGTTGGTCTGCCTGTGAAAAGAAGT
    TGATTATTTAAAAGGACAATCTGAGATTGAAATAACAACAAGGCAAATGCTTACTTATCCACCTGTATCAATccaattta
    caaatgcctgatATTTCTGACATGACCTGCTTTCTCTTTCTCCAGCTCTTCCACCTTGCTCTTCTGGAGCTGCAGTTCAG
    TCTTGGTGACGCTCAGCTCCTGTCTCTGTTCCACCATCATGTCTCTCAGCTCCTCCAGCTTAGACCAGATGTCAGGGGTC
    TTGTGGGCCGTCTCTGTAGCTTCAGCCCATGCCCCAAACAGAGAAAACAGCAGCAGAGCTACAGCACCTCTCATTCTGAA
    ACGACACATCTTCTGAAATGATCTAGTCTTTGAGGCTAGAAATTTACAGGATTAAAGTAGGAGCTATATCTGTGTGAAAT
    ATATTGTGTCTCTGTGGA

# 4 EXTRACT C1Q WITH 5’/3’ BUFFER REGIONS

Since we want to sequence the entirety of the C1q gene, we need to add
some buffer sequence outside of the 5’/3’ ends of the gene. I’ve
arbitrarily established large buffer regions:

- Left buffer: 500bp
- Right buffer: 500bp

``` bash
# Load bash variables into memory
source .bashvars

# Format region for use with pyfaidx
# Use awk to add/subtract desired buffer regions
region=$(awk '$3=="gene"' "${data_dir}/${genome_gff}" \
| grep "${sequence_ID}" \
| awk -v left_buffer="$left_buffer" -v right_buffer="$right_buffer" '{print $1":"$4 - left_buffer"-"$5 + right_buffer}')

echo "Region formatting: ${region}"
echo ""

# Extract region
${pyfaidx} "${data_dir}/${genome_fasta}" "${region}" | tee "${output_top}/${c1q_buffer_fasta}"
```

    Region formatting: NC_052339.1:19244056-19248513

    >NC_052339.1:19244056-19248513
    catgttctgttataatctccacccggcacagccaaaagaggactggccacccctcatagcctggttcctctctaggtttc
    ttcctaggttttggcctttctagggagtttttcctagccaccgtgcttctacacctgcattgcttgctgtttggggtttt
    aggctgggtttctgtacagcactttgagatatcagctgatgtaagaagggctatataaatacatttgatttgatttgata
    aggtcACTAGTTACTTATTAAATATGGATATGTCTGAACGACTAACCCTGAGACACATCAGAAAGCTGAGTTCACGTCAT
    AATATTGTGATTATAATCACAATCCCTTCTTACACCACTATAAATGTGTCTCTGTCCTTATCTGAGTCTCTGAATGTAGC
    CTATGTGTGACCTGAATGTATCTGTTACCTTAATGTAATGTAGCCTACATGTGTTACTTGAAAAATAACCCTAAATTGGT
    TCTTATTCCCAACACAGCTCTGGGATGTGTAATGGTCATATTTATTACCATACCTCTCCACTTCTTCCTCTTTCCcttta
    caatttttattttttttaaataataaagtGATGATTTATTTATTGGTTAAAGAGGAAAATCCTCACATCACAGAATGGTG
    CTGTCCTCACACTGGGAAGAGCAGGAAGCCACTGAAGGTGTTGTGGTTATCATGGCTATCATGGAGACCCTGGTCTTCAG
    GGAGACGAAGGTAGACCACATCCTCCTTCTCCAGCTCTAGTGTCAACGCATTAGATATGTACTGCCAACCCCCAACTGTA
    TTTCTTTCTACATTATATAGAACTCTTTGATTATTGTGAAACATCATGATGCCCATGAGTTGTGGGTGGCGGGATGACAT
    GGCCGTGAATCTGAAGAAGTAGACTCCTTTCACTGATGCTGTGAAGATGCCTGCATCAGAGGAGAAATAGAGTGATATCA
    GCACTGTCTGTCTCATGTTGCTGTAGCTGGAATGTTTTtaacaaataaactcagcaaaaaaagaaatgtcccttttcagg
    accctgtctttcaaagataattcgtaaaaatccaaataacttcacagatcttcattttaaagggtttaaacactgtttcc
    catgcttcttcaatgaaccataaacaattaatgaacatgcacctgtggaacggtcgttaagacactaacagcttacagac
    ggtaggcagttAGGCcaaagttatgaaaacttaggacactaaagaggcctttctactgacactgaaaaacaccaaaagaa
    agatgtccagggtccctgctcatctgcaggaacgtgccttagacatgctgcaaggaggcatgaggactacagatgtggcc
    agggcaataaattgcaatgtccgtactgtgagacacctaacacagcgctacagggagacaggacggacagctgatcgtcc
    tcgcagtgacagaccacgtgtaacaacacctgcacaggatcggtacatccgaacatcacacctgcgggacaggtaccgga
    tggcaacaacaactgcccgagttacaccaggaacgcacaatccctccatcagtgctcagactgtctgcaataggctaaga
    gaggttggactgagggcttgtaggcctgttgtcaGGCAGGTCATCACCGGCAActacgtcgcctatgggcacaaacccac
    cgtcactggaccagacagccgcggttttgtctcaccaggggtgatggtcggaattgcgtttatcgttgaaggaatgagcg
    ttacactgaggcctgtactctggagcgggatcgatttggaggtggagggtccgtcatggtctggggcggtgtatcacagc
    atcatcggactgagcttgttgtctttgcaggaaatctcaacgctgtgcgttacagggaagacatcctcctccctcatatg
    gtacccttcctgcaggctcatcctgacatgaccctccagcatgacaatgccaccagccatactgctcatgttgtgcgtga
    tttcctgcaagacaggaatgtcagtgttctgccatggccagcgaagagcccagatctcaatcccattgagcacgtctggg
    accttttggatcggagggtgagagctagggccattcgccacagaaatgtctgggaacttgcaggtgccttggtggaagag
    tggggtaacatctcacagcaagaactggcaaatctggtgcagtccatgaggaggagatgcactgcagtacttaatgcagc
    tggtggccacaccagatactgactgttacttttgattttgaccccacctttgttcagggacacattattccatttctgtt
    agtcacatgtctgtggaacttgttcagtttatgtctcagttgttgaatattgttatgttcatacaaatatttacacatgt
    taagttgactgaaaataaacgcagttgacagtgagaggacgtttcttttaaCAGTGAGAGGAGTTTAGTATAATTAAGTC
    AGGGACTAGGTAACACAAATTACCACCAAACACAATAAGTATCAAACTAAACATTGTTATAACTAAAATAAAATCAGGGT
    GACAGATTAAGTAACATAAAATATGTCATTAAAAAAAATGACCGTAgacatacagtactagtcaaaagtttggacatact
    tagaggtgtgggtttttctttatttttactattttctacattgtagcataatagtgaagacataaactatgaaataacaa
    atatggaatcatgtagtaaccaaaatagtatttaaaaaatctaaatatatttcagatttttcaaagtagccaccctttgc
    ctcgatgacagctttgcactctcttggcattctctcttcacctggaatgcttttccaagagtgtgcaaaactgtcatcaa
    ggcaaagggtgactactttaaagaatgtcaaatataaaatatattttgatttgttttatactttttgggttactacatga
    ttccatatgtgttatttaatagttttgaagtcttcactattattctacaatgtagaaaatagtaaaataaagaatccctg
    gaataagtaggtgtgtccaaacttttgactggtactgtacgtgttTGGTATCATATAGCTAGTAGTTACAAAAGTAAACT
    CATTTAAAAGCTTGAAACGTAGAATTTGTCCCTCAGAATTATGTCCACTCTGAAAATCTAGTGTCTTATCATGATTTTCA
    AAGTAAAATATAAAACAAGTTAGTGAAGTAAAATTATATACAGAGATTACTGAGAAACATGCAGCAAAATCATCCCAGAA
    GAATTACATTAAATTGCATTCATAAAAAATTATGATGAGATGACTGACTCAATATTTGGCATTGTTGGATGAATGTTATG
    CAACATTCAAATGGATTTCCCTCTTGTTAGTTAGTACCTGTAATTGGGTTGTAGGCCTTGCCGATGTTGGTGATGACTCT
    ACTGAAGATCAGTGTAGTAGCAGTACTGAATGGACCTACGTTTCCAGAGTCAGTCAAACCAGCAGAGAAGGCCACAGTTG
    GTCTGCCTGTGAAAAGAAGTTGATTATTTAAAAGGACAATCTGAGATTGAAATAACAACAAGGCAAATGCTTACTTATCC
    ACCTGTATCAATccaatttacaaatgcctgatATTTCTGACATGACCTGCTTTCTCTTTCTCCAGCTCTTCCACCTTGCT
    CTTCTGGAGCTGCAGTTCAGTCTTGGTGACGCTCAGCTCCTGTCTCTGTTCCACCATCATGTCTCTCAGCTCCTCCAGCT
    TAGACCAGATGTCAGGGGTCTTGTGGGCCGTCTCTGTAGCTTCAGCCCATGCCCCAAACAGAGAAAACAGCAGCAGAGCT
    ACAGCACCTCTCATTCTGAAACGACACATCTTCTGAAATGATCTAGTCTTTGAGGCTAGAAATTTACAGGATTAAAGTAG
    GAGCTATATCTGTGTGAAATATATTGTGTCTCTGTGGAAAAAAGAGGGTTCTTAAATATCAGGGAGTATTATTTGATTAG
    GATTACTGGTAGGATTCTTTATTTGTTTTTACTGGCAATGAGGTGAAGGTTTCTTGCGTCATTCACCAGTTGTTGTTACT
    GTATTTCGTCATACTTAATTTTGATAGTCCCCTAGATAGGGAgtgggcgtacacatcacagacaaactaaaatggtccac
    acagacagtgtggtaaagaaggaaggccaaaaagatcatcaaggacaacaaccacccgagccactgcctgttcacactgc
    tatcatccagaaggtgaggtcagtacaggtacatcaagctgggaccgagagattgagaaacagcttctatctcaaggcca
    tcagactgctaaacagcaatcattaACTGAGAGAGTCTGCTGCCCACATTGAGAACCAATCAcaggacactttaataaat
    ggatcactagtcactttaaacaatgccactttaaataatggcactttaataatgctta

# 5 PRIMER DESIGN USING [PRIMER3](https://github.com/primer3-org/primer3)

## 5.1 Design primers

Quick explanation: [Primer3](https://github.com/primer3-org/primer3)
requires a specially formatted input file. The file must be formatted
similarly to this:

    SEQUENCE_ID=${seq_id}
    SEQUENCE_TEMPLATE=${sequence}
    PRIMER_TASK=generic
    PRIMER_PICK_LEFT_PRIMER=3
    PRIMER_PICK_RIGHT_PRIMER=3
    PRIMER_OPT_SIZE=18
    PRIMER_MIN_SIZE=15
    PRIMER_MAX_SIZE=21
    PRIMER_MAX_NS_ACCEPTED=1
    PRIMER_PRODUCT_SIZE_RANGE=75-150
    P3_FILE_FLAG=1
    PRIMER_EXPLAIN_FLAG=1
    =

Values after the `=` on each line can be changed to whatever values the
user decides. The \${sequence} must be a nucletoide sequence on a single
line, with no line breaks.

The code in the chunk below uses a heredoc to write this information to
a file. Use of a heredoc allows the variables specified in the Primer3
config to expand to their actual values. Everything *between* the
following two lines gets printed (via cat) as shown and then redirected
to the indicated file (`primer3-params.txt`):

``` bash
cat << EOF > ${output_top}/primer3-params.txt
EOF
```

[Primer3](https://github.com/primer3-org/primer3) is run with the
`--format_output` to make a nice, human-readable output format.

I’ve also set [Primer3](https://github.com/primer3-org/primer3) to look
for sequencing primers and have defined the `SEQUENCE_TARGET`.

``` bash
# Load bash variables into memory
source .bashvars

# Get sequence only, by skipping the first record
# Remove newlines so sequence is on a single line
c1q_sequence=$(awk 'NR > 1' "${output_top}/${c1q_fasta}" | tr -d '\n')

# Determine length of C1q
c1q_length=${#c1q_sequence}

# Get sequence only, by skipping the first record
# Remove newlines so sequence is on a single line
c1q_buffer_sequence=$(awk 'NR > 1' "${output_top}/${c1q_buffer_fasta}" | tr -d '\n')

# Determine length of C1q with the 5'/3' buffer sequences included
c1q_buffer_length=${#c1q_buffer_sequence}

# Calculate combined length of 5'buffer region and C1q length
c1q_buffer_and_seq_length=$((left_buffer + c1q_length))


# Use heredoc to create Primer3 parameters file
cat << EOF > ${output_top}/primer3-params.txt
SEQUENCE_ID=${sequence_ID}
SEQUENCE_TEMPLATE=${c1q_buffer_sequence}
PRIMER_TASK=pick_sequencing_primers
SEQUENCE_TARGET=${left_buffer},${c1q_length}
PRIMER_MIN_TM=50
PRIMER_OPT_TM=60
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_MAX_NS_ACCEPTED=0
P3_FILE_FLAG=1
PRIMER_EXPLAIN_FLAG=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=${primer3_config}
PRIMER_NUM_RETURN=10
PRIMER_PAIR_EXPLAIN=considered 0, ok 0
=
EOF


# Run Primer3
${primer3} \
--format_output \
--output="${output_top}/primer3-primers.txt" \
"${output_top}/primer3-params.txt"

# Run Primer3 with default output for parsable results
${primer3} \
--output="${output_top}/primer3-primers-default-format.txt" \
"${output_top}/primer3-params.txt"


# Pull out any primers falling outside of target
tail -n 30 "${output_top}/primer3-primers.txt" \
| head -n -12 \
| awk -v left_buffer="$left_buffer" -v target="$c1q_buffer_and_seq_length" '$3 <= left_buffer || $3 >= target {print $0}'


### Write the left and right primer start positions to files for later use in R chunk(s).
# Pull out left primer start position
left_primer_start=$(tail -n 30 "${output_top}/primer3-primers.txt" \
| head -n -12 \
| awk -v left_buffer="$left_buffer" -v target="$c1q_buffer_and_seq_length" '$3 <= left_buffer || $3 >= target {print $0}' \
| grep "PRIMER" \
| awk 'NR == 1 {print $3}')

# Pull out right primer start position
right_primer_start=$(tail -n 30 "${output_top}/primer3-primers.txt" \
| head -n -12 \
| awk -v left_buffer="$left_buffer" -v target="$c1q_buffer_and_seq_length" '$3 <= left_buffer || $3 >= target {print $0}' \
| grep "PRIMER" \
| awk 'NR == 2 {print $3}')

echo "${left_primer_start}" > "${output_top}/left-primer-start.txt"
echo "${right_primer_start}" > "${output_top}/right-primer-start.txt"

echo ""
echo ""
echo "----------------------------------------------------------------------------------"
echo ""
echo "Primer results file:"
echo ""
echo ""
# Print the full output file
cat "${output_top}/primer3-primers.txt"

echo ""
echo ""
echo "----------------------------------------------------------------------------------"
echo ""
echo "Primer default (i.e. computer-friendly) results file:"
echo ""
echo ""

# Print the full output file
cat "${output_top}/primer3-primers-default-format.txt"
```

     1 LEFT_PRIMER        277   20   58.01   50.00    0.00   0.00    0.00 ACGACTAACCCTGAGACACA

                         start  len      tm     gc%  any_th  3'_th hairpin seq
     8 RIGHT_PRIMER      4188   21   59.22   47.62    3.32   0.00    0.00 tggccttccttctttaccaca


    ----------------------------------------------------------------------------------

    Primer results file:


    PRIMER PICKING RESULTS FOR LOC120027825

    No mispriming library specified
    Using 0-based sequence positions
    SEQUENCE SIZE: 4458
    INCLUDED REGION SIZE: 4458

    TARGETS (start, len)*: 500,3458
        0 catgttctgttataatctccacccggcacagccaaaagaggactggccacccctcatagc
                                                                      

       60 ctggttcctctctaggtttcttcctaggttttggcctttctagggagtttttcctagcca
                                                                      

      120 ccgtgcttctacacctgcattgcttgctgtttggggttttaggctgggtttctgtacagc
                                                                      

      180 actttgagatatcagctgatgtaagaagggctatataaatacatttgatttgatttgata
                                                                      

      240 aggtcACTAGTTACTTATTAAATATGGATATGTCTGAACGACTAACCCTGAGACACATCA
                                               >>>>>>>>>>>>>>>>>>>>   

      300 GAAAGCTGAGTTCACGTCATAATATTGTGATTATAATCACAATCCCTTCTTACACCACTA
                                                                      

      360 TAAATGTGTCTCTGTCCTTATCTGAGTCTCTGAATGTAGCCTATGTGTGACCTGAATGTA
                                                                      

      420 TCTGTTACCTTAATGTAATGTAGCCTACATGTGTTACTTGAAAAATAACCCTAAATTGGT
                                                                      

      480 TCTTATTCCCAACACAGCTCTGGGATGTGTAATGGTCATATTTATTACCATACCTCTCCA
                              ****************************************

      540 CTTCTTCCTCTTTCCctttacaatttttattttttttaaataataaagtGATGATTTATT
          ************************************************************

      600 TATTGGTTAAAGAGGAAAATCCTCACATCACAGAATGGTGCTGTCCTCACACTGGGAAGA
          **********************************<<<<<<<<<<<<<<<<<<<<******

      660 GCAGGAAGCCACTGAAGGTGTTGTGGTTATCATGGCTATCATGGAGACCCTGGTCTTCAG
          ************************************************************

      720 GGAGACGAAGGTAGACCACATCCTCCTTCTCCAGCTCTAGTGTCAACGCATTAGATATGT
          ************************************************************

      780 ACTGCCAACCCCCAACTGTATTTCTTTCTACATTATATAGAACTCTTTGATTATTGTGAA
          *>>>>>>>>>>>>>>>>>>>>***************************************

      840 ACATCATGATGCCCATGAGTTGTGGGTGGCGGGATGACATGGCCGTGAATCTGAAGAAGT
          ************************************************************

      900 AGACTCCTTTCACTGATGCTGTGAAGATGCCTGCATCAGAGGAGAAATAGAGTGATATCA
          ************************************************************

      960 GCACTGTCTGTCTCATGTTGCTGTAGCTGGAATGTTTTtaacaaataaactcagcaaaaa
          ************************************************************

     1020 aagaaatgtcccttttcaggaccctgtctttcaaagataattcgtaaaaatccaaataac
          ************************************************************

     1080 ttcacagatcttcattttaaagggtttaaacactgtttcccatgcttcttcaatgaacca
          ************************************************************

     1140 taaacaattaatgaacatgcacctgtggaacggtcgttaagacactaacagcttacagac
          *******************<<<<<<<<<<<<<<<<<<<<*********************

     1200 ggtaggcagttAGGCcaaagttatgaaaacttaggacactaaagaggcctttctactgac
          ************************************************************

     1260 actgaaaaacaccaaaagaaagatgtccagggtccctgctcatctgcaggaacgtgcctt
          *******************************************>>>>>>>>>>>>>>>>>

     1320 agacatgctgcaaggaggcatgaggactacagatgtggccagggcaataaattgcaatgt
          >>>*********************************************************

     1380 ccgtactgtgagacacctaacacagcgctacagggagacaggacggacagctgatcgtcc
          ************************************************************

     1440 tcgcagtgacagaccacgtgtaacaacacctgcacaggatcggtacatccgaacatcaca
          ************************************************************

     1500 cctgcgggacaggtaccggatggcaacaacaactgcccgagttacaccaggaacgcacaa
          ************************************************************

     1560 tccctccatcagtgctcagactgtctgcaataggctaagagaggttggactgagggcttg
          ************************************************************

     1620 taggcctgttgtcaGGCAGGTCATCACCGGCAActacgtcgcctatgggcacaaacccac
          ****************<<<<<<<<<<<<<<<<<<<<************************

     1680 cgtcactggaccagacagccgcggttttgtctcaccaggggtgatggtcggaattgcgtt
          ************************************************************

     1740 tatcgttgaaggaatgagcgttacactgaggcctgtactctggagcgggatcgatttgga
          **************************************>>>>>>>>>>>>>>>>>>>>**

     1800 ggtggagggtccgtcatggtctggggcggtgtatcacagcatcatcggactgagcttgtt
          ************************************************************

     1860 gtctttgcaggaaatctcaacgctgtgcgttacagggaagacatcctcctccctcatatg
          ************************************************************

     1920 gtacccttcctgcaggctcatcctgacatgaccctccagcatgacaatgccaccagccat
          ************************************************************

     1980 actgctcatgttgtgcgtgatttcctgcaagacaggaatgtcagtgttctgccatggcca
          ************************************************************

     2040 gcgaagagcccagatctcaatcccattgagcacgtctgggaccttttggatcggagggtg
          ************************************************************

     2100 agagctagggccattcgccacagaaatgtctgggaacttgcaggtgccttggtggaagag
          ************************************************************

     2160 tggggtaacatctcacagcaagaactggcaaatctggtgcagtccatgaggaggagatgc
          **********<<<<<<<<<<<<<<<<<<<<******************************

     2220 actgcagtacttaatgcagctggtggccacaccagatactgactgttacttttgattttg
          ************************************************************

     2280 accccacctttgttcagggacacattattccatttctgttagtcacatgtctgtggaact
          **>>>>>>>>>>>>>>>>>>>>**************************************

     2340 tgttcagtttatgtctcagttgttgaatattgttatgttcatacaaatatttacacatgt
          ************************************************************

     2400 taagttgactgaaaataaacgcagttgacagtgagaggacgtttcttttaaCAGTGAGAG
          ************************************************************

     2460 GAGTTTAGTATAATTAAGTCAGGGACTAGGTAACACAAATTACCACCAAACACAATAAGT
          ************************************************************

     2520 ATCAAACTAAACATTGTTATAACTAAAATAAAATCAGGGTGACAGATTAAGTAACATAAA
          ************************************************************

     2580 ATATGTCATTAAAAAAAATGACCGTAgacatacagtactagtcaaaagtttggacatact
          *********************************************************<<<

     2640 tagaggtgtgggtttttctttatttttactattttctacattgtagcataatagtgaaga
          <<<<<<<<<<<<<<<<<<<*****************************************

     2700 cataaactatgaaataacaaatatggaatcatgtagtaaccaaaatagtatttaaaaaat
          ************************************************************

     2760 ctaaatatatttcagatttttcaaagtagccaccctttgcctcgatgacagctttgcact
          *******************************>>>>>>>>>>>>>>>>>>>>*********

     2820 ctcttggcattctctcttcacctggaatgcttttccaagagtgtgcaaaactgtcatcaa
          ************************************************************

     2880 ggcaaagggtgactactttaaagaatgtcaaatataaaatatattttgatttgttttata
          ************************************************************

     2940 ctttttgggttactacatgattccatatgtgttatttaatagttttgaagtcttcactat
          ************************************************************

     3000 tattctacaatgtagaaaatagtaaaataaagaatccctggaataagtaggtgtgtccaa
          ************************************************************

     3060 acttttgactggtactgtacgtgttTGGTATCATATAGCTAGTAGTTACAAAAGTAAACT
          ************************************************************

     3120 CATTTAAAAGCTTGAAACGTAGAATTTGTCCCTCAGAATTATGTCCACTCTGAAAATCTA
          **************************<<<<<<<<<<<<<<<<<<<<<<<***********

     3180 GTGTCTTATCATGATTTTCAAAGTAAAATATAAAACAAGTTAGTGAAGTAAAATTATATA
          ************************************************************

     3240 CAGAGATTACTGAGAAACATGCAGCAAAATCATCCCAGAAGAATTACATTAAATTGCATT
          ********************>>>>>>>>>>>>>>>>>>>>>>>>>***************

     3300 CATAAAAAATTATGATGAGATGACTGACTCAATATTTGGCATTGTTGGATGAATGTTATG
          ************************************************************

     3360 CAACATTCAAATGGATTTCCCTCTTGTTAGTTAGTACCTGTAATTGGGTTGTAGGCCTTG
          ************************************************************

     3420 CCGATGTTGGTGATGACTCTACTGAAGATCAGTGTAGTAGCAGTACTGAATGGACCTACG
          ************************************************************

     3480 TTTCCAGAGTCAGTCAAACCAGCAGAGAAGGCCACAGTTGGTCTGCCTGTGAAAAGAAGT
          ************************************************************

     3540 TGATTATTTAAAAGGACAATCTGAGATTGAAATAACAACAAGGCAAATGCTTACTTATCC
          ************************************************************

     3600 ACCTGTATCAATccaatttacaaatgcctgatATTTCTGACATGACCTGCTTTCTCTTTC
          ************************************************************

     3660 TCCAGCTCTTCCACCTTGCTCTTCTGGAGCTGCAGTTCAGTCTTGGTGACGCTCAGCTCC
          **<<<<<<<<<<<<<<<<<<<<**************************************

     3720 TGTCTCTGTTCCACCATCATGTCTCTCAGCTCCTCCAGCTTAGACCAGATGTCAGGGGTC
          ***********************************************>>>>>>>>>>>>>

     3780 TTGTGGGCCGTCTCTGTAGCTTCAGCCCATGCCCCAAACAGAGAAAACAGCAGCAGAGCT
          >>>>>>>*****************************************************

     3840 ACAGCACCTCTCATTCTGAAACGACACATCTTCTGAAATGATCTAGTCTTTGAGGCTAGA
          ************************************************************

     3900 AATTTACAGGATTAAAGTAGGAGCTATATCTGTGTGAAATATATTGTGTCTCTGTGGAAA
          **********************************************************  

     3960 AAAGAGGGTTCTTAAATATCAGGGAGTATTATTTGATTAGGATTACTGGTAGGATTCTTT
                                                                      

     4020 ATTTGTTTTTACTGGCAATGAGGTGAAGGTTTCTTGCGTCATTCACCAGTTGTTGTTACT
                                                                      

     4080 GTATTTCGTCATACTTAATTTTGATAGTCCCCTAGATAGGGAgtgggcgtacacatcaca
                                                                      

     4140 gacaaactaaaatggtccacacagacagtgtggtaaagaaggaaggccaaaaagatcatc
                                      <<<<<<<<<<<<<<<<<<<<<           

     4200 aaggacaacaaccacccgagccactgcctgttcacactgctatcatccagaaggtgaggt
                                                                      

     4260 cagtacaggtacatcaagctgggaccgagagattgagaaacagcttctatctcaaggcca
                                                                      

     4320 tcagactgctaaacagcaatcattaACTGAGAGAGTCTGCTGCCCACATTGAGAACCAAT
                                                                      

     4380 CAcaggacactttaataaatggatcactagtcactttaaacaatgccactttaaataatg
                                                                      

     4440 gcactttaataatgctta
                            

    KEYS (in order of precedence):
    ****** target
    >>>>>> left primer
    <<<<<< right primer
    ^^^^^^ left primer / right primer overlap

                        start  len      tm     gc%  any_th  3'_th hairpin seq
     1 LEFT_PRIMER        277   20   58.01   50.00    0.00   0.00    0.00 ACGACTAACCCTGAGACACA
     2 LEFT_PRIMER        781   20   59.96   55.00    0.00   0.00    0.00 CTGCCAACCCCCAACTGTAT
     3 LEFT_PRIMER       1303   20   60.04   55.00   14.25   0.00   43.28 ctgcaggaacgtgccttaga
     4 LEFT_PRIMER       1778   20   59.90   55.00    0.00   0.00   39.00 tctggagcgggatcgatttg
     5 LEFT_PRIMER       2282   20   59.82   55.00    0.00   0.00   42.25 cccacctttgttcagggaca
     6 LEFT_PRIMER       2791   20   60.04   55.00    0.00   0.00    0.00 accctttgcctcgatgacag
     7 LEFT_PRIMER       3260   25   60.63   40.00    0.00   0.00    0.00 GCAGCAAAATCATCCCAGAAGAATT
     8 LEFT_PRIMER       3767   20   60.04   60.00    0.00   0.00    0.00 GATGTCAGGGGTCTTGTGGG

                         start  len      tm     gc%  any_th  3'_th hairpin seq
     1 RIGHT_PRIMER       653   20   59.96   55.00    0.00   0.00   35.43 CAGTGTGAGGACAGCACCAT
     2 RIGHT_PRIMER      1178   20   59.97   55.00   11.27   0.00   41.51 taacgaccgttccacaggtg
     3 RIGHT_PRIMER      1655   20   59.75   55.00    1.90   0.00   40.58 tagTTGCCGGTGATGACCTG
     4 RIGHT_PRIMER      2189   20   59.53   50.00    0.00   0.00   34.84 tgccagttcttgctgtgaga
     5 RIGHT_PRIMER      2658   22   57.21   40.91    0.00   0.00    0.00 agaaaaacccacacctctaagt
     6 RIGHT_PRIMER      3168   23   58.90   43.48    0.00   0.00    0.00 AGTGGACATAATTCTGAGGGACA
     7 RIGHT_PRIMER      3681   20   59.68   55.00    0.00   0.00   35.36 AGAGCAAGGTGGAAGAGCTG
     8 RIGHT_PRIMER      4188   21   59.22   47.62    3.32   0.00    0.00 tggccttccttctttaccaca

    Statistics
             con   too    in    in   not          no    tm    tm   high  high  high        high      
             sid  many   tar  excl    ok   bad    GC   too   too any_th 3'_th hair-  poly   end      
            ered    Ns   get   reg   reg   GC% clamp   low  high  compl compl   pin     X  stab    ok
    Left    2518     0     0     0     0   162     0   221   612      0     0   176    10     0     8
    Right   2560     0     0     0     0   108     0   205   530      0     0   221     0     0     8
    Pair Stats:
    considered 0, ok 0
    libprimer3 release 2.6.1




    ----------------------------------------------------------------------------------

    Primer default (i.e. computer-friendly) results file:


    SEQUENCE_ID=LOC120027825
    SEQUENCE_TEMPLATE=catgttctgttataatctccacccggcacagccaaaagaggactggccacccctcatagcctggttcctctctaggtttcttcctaggttttggcctttctagggagtttttcctagccaccgtgcttctacacctgcattgcttgctgtttggggttttaggctgggtttctgtacagcactttgagatatcagctgatgtaagaagggctatataaatacatttgatttgatttgataaggtcACTAGTTACTTATTAAATATGGATATGTCTGAACGACTAACCCTGAGACACATCAGAAAGCTGAGTTCACGTCATAATATTGTGATTATAATCACAATCCCTTCTTACACCACTATAAATGTGTCTCTGTCCTTATCTGAGTCTCTGAATGTAGCCTATGTGTGACCTGAATGTATCTGTTACCTTAATGTAATGTAGCCTACATGTGTTACTTGAAAAATAACCCTAAATTGGTTCTTATTCCCAACACAGCTCTGGGATGTGTAATGGTCATATTTATTACCATACCTCTCCACTTCTTCCTCTTTCCctttacaatttttattttttttaaataataaagtGATGATTTATTTATTGGTTAAAGAGGAAAATCCTCACATCACAGAATGGTGCTGTCCTCACACTGGGAAGAGCAGGAAGCCACTGAAGGTGTTGTGGTTATCATGGCTATCATGGAGACCCTGGTCTTCAGGGAGACGAAGGTAGACCACATCCTCCTTCTCCAGCTCTAGTGTCAACGCATTAGATATGTACTGCCAACCCCCAACTGTATTTCTTTCTACATTATATAGAACTCTTTGATTATTGTGAAACATCATGATGCCCATGAGTTGTGGGTGGCGGGATGACATGGCCGTGAATCTGAAGAAGTAGACTCCTTTCACTGATGCTGTGAAGATGCCTGCATCAGAGGAGAAATAGAGTGATATCAGCACTGTCTGTCTCATGTTGCTGTAGCTGGAATGTTTTtaacaaataaactcagcaaaaaaagaaatgtcccttttcaggaccctgtctttcaaagataattcgtaaaaatccaaataacttcacagatcttcattttaaagggtttaaacactgtttcccatgcttcttcaatgaaccataaacaattaatgaacatgcacctgtggaacggtcgttaagacactaacagcttacagacggtaggcagttAGGCcaaagttatgaaaacttaggacactaaagaggcctttctactgacactgaaaaacaccaaaagaaagatgtccagggtccctgctcatctgcaggaacgtgccttagacatgctgcaaggaggcatgaggactacagatgtggccagggcaataaattgcaatgtccgtactgtgagacacctaacacagcgctacagggagacaggacggacagctgatcgtcctcgcagtgacagaccacgtgtaacaacacctgcacaggatcggtacatccgaacatcacacctgcgggacaggtaccggatggcaacaacaactgcccgagttacaccaggaacgcacaatccctccatcagtgctcagactgtctgcaataggctaagagaggttggactgagggcttgtaggcctgttgtcaGGCAGGTCATCACCGGCAActacgtcgcctatgggcacaaacccaccgtcactggaccagacagccgcggttttgtctcaccaggggtgatggtcggaattgcgtttatcgttgaaggaatgagcgttacactgaggcctgtactctggagcgggatcgatttggaggtggagggtccgtcatggtctggggcggtgtatcacagcatcatcggactgagcttgttgtctttgcaggaaatctcaacgctgtgcgttacagggaagacatcctcctccctcatatggtacccttcctgcaggctcatcctgacatgaccctccagcatgacaatgccaccagccatactgctcatgttgtgcgtgatttcctgcaagacaggaatgtcagtgttctgccatggccagcgaagagcccagatctcaatcccattgagcacgtctgggaccttttggatcggagggtgagagctagggccattcgccacagaaatgtctgggaacttgcaggtgccttggtggaagagtggggtaacatctcacagcaagaactggcaaatctggtgcagtccatgaggaggagatgcactgcagtacttaatgcagctggtggccacaccagatactgactgttacttttgattttgaccccacctttgttcagggacacattattccatttctgttagtcacatgtctgtggaacttgttcagtttatgtctcagttgttgaatattgttatgttcatacaaatatttacacatgttaagttgactgaaaataaacgcagttgacagtgagaggacgtttcttttaaCAGTGAGAGGAGTTTAGTATAATTAAGTCAGGGACTAGGTAACACAAATTACCACCAAACACAATAAGTATCAAACTAAACATTGTTATAACTAAAATAAAATCAGGGTGACAGATTAAGTAACATAAAATATGTCATTAAAAAAAATGACCGTAgacatacagtactagtcaaaagtttggacatacttagaggtgtgggtttttctttatttttactattttctacattgtagcataatagtgaagacataaactatgaaataacaaatatggaatcatgtagtaaccaaaatagtatttaaaaaatctaaatatatttcagatttttcaaagtagccaccctttgcctcgatgacagctttgcactctcttggcattctctcttcacctggaatgcttttccaagagtgtgcaaaactgtcatcaaggcaaagggtgactactttaaagaatgtcaaatataaaatatattttgatttgttttatactttttgggttactacatgattccatatgtgttatttaatagttttgaagtcttcactattattctacaatgtagaaaatagtaaaataaagaatccctggaataagtaggtgtgtccaaacttttgactggtactgtacgtgttTGGTATCATATAGCTAGTAGTTACAAAAGTAAACTCATTTAAAAGCTTGAAACGTAGAATTTGTCCCTCAGAATTATGTCCACTCTGAAAATCTAGTGTCTTATCATGATTTTCAAAGTAAAATATAAAACAAGTTAGTGAAGTAAAATTATATACAGAGATTACTGAGAAACATGCAGCAAAATCATCCCAGAAGAATTACATTAAATTGCATTCATAAAAAATTATGATGAGATGACTGACTCAATATTTGGCATTGTTGGATGAATGTTATGCAACATTCAAATGGATTTCCCTCTTGTTAGTTAGTACCTGTAATTGGGTTGTAGGCCTTGCCGATGTTGGTGATGACTCTACTGAAGATCAGTGTAGTAGCAGTACTGAATGGACCTACGTTTCCAGAGTCAGTCAAACCAGCAGAGAAGGCCACAGTTGGTCTGCCTGTGAAAAGAAGTTGATTATTTAAAAGGACAATCTGAGATTGAAATAACAACAAGGCAAATGCTTACTTATCCACCTGTATCAATccaatttacaaatgcctgatATTTCTGACATGACCTGCTTTCTCTTTCTCCAGCTCTTCCACCTTGCTCTTCTGGAGCTGCAGTTCAGTCTTGGTGACGCTCAGCTCCTGTCTCTGTTCCACCATCATGTCTCTCAGCTCCTCCAGCTTAGACCAGATGTCAGGGGTCTTGTGGGCCGTCTCTGTAGCTTCAGCCCATGCCCCAAACAGAGAAAACAGCAGCAGAGCTACAGCACCTCTCATTCTGAAACGACACATCTTCTGAAATGATCTAGTCTTTGAGGCTAGAAATTTACAGGATTAAAGTAGGAGCTATATCTGTGTGAAATATATTGTGTCTCTGTGGAAAAAAGAGGGTTCTTAAATATCAGGGAGTATTATTTGATTAGGATTACTGGTAGGATTCTTTATTTGTTTTTACTGGCAATGAGGTGAAGGTTTCTTGCGTCATTCACCAGTTGTTGTTACTGTATTTCGTCATACTTAATTTTGATAGTCCCCTAGATAGGGAgtgggcgtacacatcacagacaaactaaaatggtccacacagacagtgtggtaaagaaggaaggccaaaaagatcatcaaggacaacaaccacccgagccactgcctgttcacactgctatcatccagaaggtgaggtcagtacaggtacatcaagctgggaccgagagattgagaaacagcttctatctcaaggccatcagactgctaaacagcaatcattaACTGAGAGAGTCTGCTGCCCACATTGAGAACCAATCAcaggacactttaataaatggatcactagtcactttaaacaatgccactttaaataatggcactttaataatgctta
    PRIMER_TASK=pick_sequencing_primers
    SEQUENCE_TARGET=500,3458
    PRIMER_MIN_TM=50
    PRIMER_OPT_TM=60
    PRIMER_PICK_LEFT_PRIMER=1
    PRIMER_PICK_RIGHT_PRIMER=1
    PRIMER_OPT_SIZE=20
    PRIMER_MIN_SIZE=18
    PRIMER_MAX_SIZE=25
    PRIMER_MAX_NS_ACCEPTED=0
    P3_FILE_FLAG=1
    PRIMER_EXPLAIN_FLAG=1
    PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/home/shared/primer3-2.6.1/src/primer3_config
    PRIMER_NUM_RETURN=10
    PRIMER_PAIR_EXPLAIN=considered 0, ok 0
    PRIMER_LEFT_EXPLAIN=sequencing locations 8, considered 2518, GC content failed 162, low tm 221, high tm 612, high hairpin stability 176, long poly-x seq 10, ok 8
    PRIMER_RIGHT_EXPLAIN=sequencing locations 8, considered 2560, GC content failed 108, low tm 205, high tm 530, high hairpin stability 221, ok 8
    PRIMER_PAIR_EXPLAIN=considered 0, ok 0
    PRIMER_LEFT_NUM_RETURNED=8
    PRIMER_RIGHT_NUM_RETURNED=8
    PRIMER_INTERNAL_NUM_RETURNED=0
    PRIMER_PAIR_NUM_RETURNED=0
    PRIMER_LEFT_0_PENALTY=1.991121
    PRIMER_RIGHT_0_PENALTY=0.035881
    PRIMER_LEFT_0_SEQUENCE=ACGACTAACCCTGAGACACA
    PRIMER_RIGHT_0_SEQUENCE=CAGTGTGAGGACAGCACCAT
    PRIMER_LEFT_0=277,20
    PRIMER_RIGHT_0=653,20
    PRIMER_LEFT_0_TM=58.009
    PRIMER_RIGHT_0_TM=59.964
    PRIMER_LEFT_0_GC_PERCENT=50.000
    PRIMER_RIGHT_0_GC_PERCENT=55.000
    PRIMER_LEFT_0_SELF_ANY_TH=0.00
    PRIMER_RIGHT_0_SELF_ANY_TH=0.00
    PRIMER_LEFT_0_SELF_END_TH=0.00
    PRIMER_RIGHT_0_SELF_END_TH=0.00
    PRIMER_LEFT_0_HAIRPIN_TH=0.00
    PRIMER_RIGHT_0_HAIRPIN_TH=35.43
    PRIMER_LEFT_0_END_STABILITY=3.7200
    PRIMER_RIGHT_0_END_STABILITY=3.5500
    PRIMER_LEFT_1_PENALTY=0.039618
    PRIMER_RIGHT_1_PENALTY=0.032403
    PRIMER_LEFT_1_SEQUENCE=CTGCCAACCCCCAACTGTAT
    PRIMER_RIGHT_1_SEQUENCE=taacgaccgttccacaggtg
    PRIMER_LEFT_1=781,20
    PRIMER_RIGHT_1=1178,20
    PRIMER_LEFT_1_TM=59.960
    PRIMER_RIGHT_1_TM=59.968
    PRIMER_LEFT_1_GC_PERCENT=55.000
    PRIMER_RIGHT_1_GC_PERCENT=55.000
    PRIMER_LEFT_1_SELF_ANY_TH=0.00
    PRIMER_RIGHT_1_SELF_ANY_TH=11.27
    PRIMER_LEFT_1_SELF_END_TH=0.00
    PRIMER_RIGHT_1_SELF_END_TH=0.00
    PRIMER_LEFT_1_HAIRPIN_TH=0.00
    PRIMER_RIGHT_1_HAIRPIN_TH=41.51
    PRIMER_LEFT_1_END_STABILITY=2.2900
    PRIMER_RIGHT_1_END_STABILITY=4.0000
    PRIMER_LEFT_2_PENALTY=0.037326
    PRIMER_RIGHT_2_PENALTY=0.249120
    PRIMER_LEFT_2_SEQUENCE=ctgcaggaacgtgccttaga
    PRIMER_RIGHT_2_SEQUENCE=tagTTGCCGGTGATGACCTG
    PRIMER_LEFT_2=1303,20
    PRIMER_RIGHT_2=1655,20
    PRIMER_LEFT_2_TM=60.037
    PRIMER_RIGHT_2_TM=59.751
    PRIMER_LEFT_2_GC_PERCENT=55.000
    PRIMER_RIGHT_2_GC_PERCENT=55.000
    PRIMER_LEFT_2_SELF_ANY_TH=14.25
    PRIMER_RIGHT_2_SELF_ANY_TH=1.90
    PRIMER_LEFT_2_SELF_END_TH=0.00
    PRIMER_RIGHT_2_SELF_END_TH=0.00
    PRIMER_LEFT_2_HAIRPIN_TH=43.28
    PRIMER_RIGHT_2_HAIRPIN_TH=40.58
    PRIMER_LEFT_2_END_STABILITY=2.1000
    PRIMER_RIGHT_2_END_STABILITY=4.0000
    PRIMER_LEFT_3_PENALTY=0.104855
    PRIMER_RIGHT_3_PENALTY=0.469494
    PRIMER_LEFT_3_SEQUENCE=tctggagcgggatcgatttg
    PRIMER_RIGHT_3_SEQUENCE=tgccagttcttgctgtgaga
    PRIMER_LEFT_3=1778,20
    PRIMER_RIGHT_3=2189,20
    PRIMER_LEFT_3_TM=59.895
    PRIMER_RIGHT_3_TM=59.531
    PRIMER_LEFT_3_GC_PERCENT=55.000
    PRIMER_RIGHT_3_GC_PERCENT=50.000
    PRIMER_LEFT_3_SELF_ANY_TH=0.00
    PRIMER_RIGHT_3_SELF_ANY_TH=0.00
    PRIMER_LEFT_3_SELF_END_TH=0.00
    PRIMER_RIGHT_3_SELF_END_TH=0.00
    PRIMER_LEFT_3_HAIRPIN_TH=39.00
    PRIMER_RIGHT_3_HAIRPIN_TH=34.84
    PRIMER_LEFT_3_END_STABILITY=2.3200
    PRIMER_RIGHT_3_END_STABILITY=3.2700
    PRIMER_LEFT_4_PENALTY=0.184933
    PRIMER_RIGHT_4_PENALTY=4.790718
    PRIMER_LEFT_4_SEQUENCE=cccacctttgttcagggaca
    PRIMER_RIGHT_4_SEQUENCE=agaaaaacccacacctctaagt
    PRIMER_LEFT_4=2282,20
    PRIMER_RIGHT_4=2658,22
    PRIMER_LEFT_4_TM=59.815
    PRIMER_RIGHT_4_TM=57.209
    PRIMER_LEFT_4_GC_PERCENT=55.000
    PRIMER_RIGHT_4_GC_PERCENT=40.909
    PRIMER_LEFT_4_SELF_ANY_TH=0.00
    PRIMER_RIGHT_4_SELF_ANY_TH=0.00
    PRIMER_LEFT_4_SELF_END_TH=0.00
    PRIMER_RIGHT_4_SELF_END_TH=0.00
    PRIMER_LEFT_4_HAIRPIN_TH=42.25
    PRIMER_RIGHT_4_HAIRPIN_TH=0.00
    PRIMER_LEFT_4_END_STABILITY=4.0200
    PRIMER_RIGHT_4_END_STABILITY=2.2400
    PRIMER_LEFT_5_PENALTY=0.036010
    PRIMER_RIGHT_5_PENALTY=4.104086
    PRIMER_LEFT_5_SEQUENCE=accctttgcctcgatgacag
    PRIMER_RIGHT_5_SEQUENCE=AGTGGACATAATTCTGAGGGACA
    PRIMER_LEFT_5=2791,20
    PRIMER_RIGHT_5=3168,23
    PRIMER_LEFT_5_TM=60.036
    PRIMER_RIGHT_5_TM=58.896
    PRIMER_LEFT_5_GC_PERCENT=55.000
    PRIMER_RIGHT_5_GC_PERCENT=43.478
    PRIMER_LEFT_5_SELF_ANY_TH=0.00
    PRIMER_RIGHT_5_SELF_ANY_TH=0.00
    PRIMER_LEFT_5_SELF_END_TH=0.00
    PRIMER_RIGHT_5_SELF_END_TH=0.00
    PRIMER_LEFT_5_HAIRPIN_TH=0.00
    PRIMER_RIGHT_5_HAIRPIN_TH=0.00
    PRIMER_LEFT_5_END_STABILITY=3.5100
    PRIMER_RIGHT_5_END_STABILITY=4.0200
    PRIMER_LEFT_6_PENALTY=5.626116
    PRIMER_RIGHT_6_PENALTY=0.324285
    PRIMER_LEFT_6_SEQUENCE=GCAGCAAAATCATCCCAGAAGAATT
    PRIMER_RIGHT_6_SEQUENCE=AGAGCAAGGTGGAAGAGCTG
    PRIMER_LEFT_6=3260,25
    PRIMER_RIGHT_6=3681,20
    PRIMER_LEFT_6_TM=60.626
    PRIMER_RIGHT_6_TM=59.676
    PRIMER_LEFT_6_GC_PERCENT=40.000
    PRIMER_RIGHT_6_GC_PERCENT=55.000
    PRIMER_LEFT_6_SELF_ANY_TH=0.00
    PRIMER_RIGHT_6_SELF_ANY_TH=0.00
    PRIMER_LEFT_6_SELF_END_TH=0.00
    PRIMER_RIGHT_6_SELF_END_TH=0.00
    PRIMER_LEFT_6_HAIRPIN_TH=0.00
    PRIMER_RIGHT_6_HAIRPIN_TH=35.36
    PRIMER_LEFT_6_END_STABILITY=2.1700
    PRIMER_RIGHT_6_END_STABILITY=4.2400
    PRIMER_LEFT_7_PENALTY=0.035056
    PRIMER_RIGHT_7_PENALTY=1.777878
    PRIMER_LEFT_7_SEQUENCE=GATGTCAGGGGTCTTGTGGG
    PRIMER_RIGHT_7_SEQUENCE=tggccttccttctttaccaca
    PRIMER_LEFT_7=3767,20
    PRIMER_RIGHT_7=4188,21
    PRIMER_LEFT_7_TM=60.035
    PRIMER_RIGHT_7_TM=59.222
    PRIMER_LEFT_7_GC_PERCENT=60.000
    PRIMER_RIGHT_7_GC_PERCENT=47.619
    PRIMER_LEFT_7_SELF_ANY_TH=0.00
    PRIMER_RIGHT_7_SELF_ANY_TH=3.32
    PRIMER_LEFT_7_SELF_END_TH=0.00
    PRIMER_RIGHT_7_SELF_END_TH=0.00
    PRIMER_LEFT_7_HAIRPIN_TH=0.00
    PRIMER_RIGHT_7_HAIRPIN_TH=0.00
    PRIMER_LEFT_7_END_STABILITY=4.6100
    PRIMER_RIGHT_7_END_STABILITY=4.1700
    =

Looks like we have one primer pair where both are outside of the C1q
gene sequence.

## 5.2 Amplicon length

The primer pair will produced an amplicon length of 3911bp.

Next, we need to test them against the entire genome to assess
specificity.

# 6 SPLIT GENOME

For some reason, I feel like the EMBOSS PrimerSearch tool will *not*
work on mulit-FastA files (or, it makes the results more difficult to
decipher?), but I’m just going off of memory and repeating a process
I’ve done previously.

This will split the genome multi-FastA file into separate FastA files.

## 6.1 Split mulit-FastA file in to individual FastA files with PyFaidx

``` bash
# Load bash variables into memory
source .bashvars

# Make directory if it doesn't exist
mkdir --parents ${genome_fasta_splits_dir}

cd ${genome_fasta_splits_dir}

# Count sequences in FastA
echo "-------------------------------------------------------------------"
echo "NUMBER OF SEQUENCES IN ORIGINAL FASTA"
grep -c "^>" ../${genome_fasta}
echo "-------------------------------------------------------------------"
echo ""
echo ""



# Split FastA
${pyfaidx} \
--split-files \
../${genome_fasta}

# Count number of individual FastA files
echo "-------------------------------------------------------------------"
echo "NUMBER OF INDIVIDUAL FASTA FILES"
ls -1 | wc -l
echo "-------------------------------------------------------------------"
```

    -------------------------------------------------------------------
    NUMBER OF SEQUENCES IN ORIGINAL FASTA
    4121
    -------------------------------------------------------------------


    -------------------------------------------------------------------
    NUMBER OF INDIVIDUAL FASTA FILES
    4122
    -------------------------------------------------------------------

# 7 PRIMER SEARCH WITH [EMBOSS PRIMERSEARCH](https://emboss.sourceforge.net/apps/cvs/emboss/apps/primersearch.html)

This will run [EMBOSS
PrimerSearch](https://emboss.sourceforge.net/apps/cvs/emboss/apps/primersearch.html)
against the genome to assess primer specificity.

## 7.1 Create EMBOSS PrimerSearch Primers File

Create a tab-delimited file to use with EMBOSS PrimerSearch.

``` bash
# Load bash variables into memory
source .bashvars

# Get primer info from Primer3 default format output
seq_id=$(grep "SEQUENCE_ID=" "${output_top}/primer3-primers-default-format.txt" | sed 's/SEQUENCE_ID=//')

left_primer=$(grep "PRIMER_LEFT_0_SEQUENCE=" "${output_top}/primer3-primers-default-format.txt" | sed 's/PRIMER_LEFT_0_SEQUENCE=//' | tr '[:lower:]' '[:upper:]')

right_primer=$(grep "PRIMER_RIGHT_7_SEQUENCE=" "${output_top}/primer3-primers-default-format.txt" | sed 's/PRIMER_RIGHT_7_SEQUENCE=//' | tr '[:lower:]' '[:upper:]')

# Create EMBOSS primer file
printf "%s\t%s\t%s\t\n" "${seq_id}" "${left_primer}" "${right_primer}" | tee "${output_top}/emboss-primers.txt"
```

    LOC120027825    ACGACTAACCCTGAGACACA    TGGCCTTCCTTCTTTACCACA   

## 7.2 Run EMBOSS PrimerSearch

This will run EMBOSS PrimerSearch and allow for a 10% mismatch in primer
annealing sites (`${primersearch} -auto ${fasta} ${primers} 10`).

Afterwards, the resulting output files (`*.primersearch`) will be
searched for the term `Amplimer`, indicating a PCR product would be
produced. If no amplimer is identified in a `*.primersearch` file, then
that file is deleted. This should leave just the results in which primer
matches were identified.

``` bash
# Load bash variables into memory
source .bashvars

cd ${genome_fasta_splits_dir}

primers="../../../../output/emboss-primers.txt"


time \
for fasta in *.fna
  do
  # Remove path from FastA filename
  fasta_no_path=$(echo ${fasta##*/})
  
  # Remove file extension from FastA filename
  fasta_no_ext=$(echo ${fasta_no_path%%.*})
  
  # Convert filename to lowercase
  # Will be used for output from EMBOSS PrimerSearch
  fasta_no_ext_lower=$(echo ${fasta_no_ext} | tr '[:upper:]' '[:lower:]')
  
  ###### Run EMBOSS PrimerSearch on all FastA files ########
  # Allows for a 10% mismatch
  ${primersearch} -auto ${fasta} ${primers} 10
  ##### END EMBOSS ##########
  
  # Find EMBOSS PrimerSearch output files with primer matches
  # Remove those without a match
  if ! grep --quiet "Amplimer" "${fasta_no_ext_lower}.primersearch"
    then rm ${fasta_no_ext_lower}.primersearch
  fi
done
```

    real    5m28.518s
    user    4m14.155s
    sys 1m45.833s

## 7.3 Check primer matches

This will print the contents of any remaining `*.primersearch` output
files (i.e. sequences with primer matches)

``` bash
# Load bash variables into memory
source .bashvars

cd ${genome_fasta_splits_dir}

# Check contents of files with matches
for file in *.primersearch
  do
  echo "FILE: ${file}"
  echo ""
  cat ${file}
  echo "----------------------------------"
  echo ""
done
```

    FILE: nc_052339.primersearch


    Primer name LOC120027825
    Amplimer 1
        Sequence: NC_052339.1  
        
        ACGACTAACCCTGAGACACA hits forward strand at 19244333 with 0 mismatches
        TGGCCTTCCTTCTTTACCACA hits reverse strand at [18091675] with 0 mismatches
        Amplimer length: 3912 bp
    ----------------------------------

# 8 SUMMARY

EMBOSS PrimerSearch has identified only a single location in the genome
where the forward and reverse primers match, as well as in the expected
location on chromosome `NC_052339`. This means our primers are specific
to our desired target in the genome.

# 9 CITATIONS

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-koressaar2007" class="csl-entry">

Koressaar, Triinu, and Maido Remm. 2007. “Enhancements and Modifications
of Primer Design Program Primer3.” *Bioinformatics* 23 (10): 1289–91.
<https://doi.org/10.1093/bioinformatics/btm091>.

</div>

<div id="ref-rice2000" class="csl-entry">

Rice, Peter, Ian Longden, and Alan Bleasby. 2000. “EMBOSS: The European
Molecular Biology Open Software Suite.” *Trends in Genetics* 16 (6):
276–77. <https://doi.org/10.1016/s0168-9525(00)02024-2>.

</div>

<div id="ref-shirley2015" class="csl-entry">

Shirley, Matthew D, Zhaorong Ma, Brent S Pedersen, and Sarah J Wheelan.
2015. “Efficient "Pythonic" Access to FASTA Files Using Pyfaidx.”
<http://dx.doi.org/10.7287/peerj.preprints.970v1>.

</div>

<div id="ref-untergasser2012" class="csl-entry">

Untergasser, Andreas, Ioana Cutcutache, Triinu Koressaar, Jian Ye, Brant
C. Faircloth, Maido Remm, and Steven G. Rozen. 2012. “Primer3new
Capabilities and Interfaces.” *Nucleic Acids Research* 40 (15): e115–15.
<https://doi.org/10.1093/nar/gks596>.

</div>

</div>
