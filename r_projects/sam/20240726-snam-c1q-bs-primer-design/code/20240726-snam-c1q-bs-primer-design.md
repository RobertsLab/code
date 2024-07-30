20240726-snam-c1q-bs-primer-design
================
Sam White
2024-07-25

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
  - <a href="#51-bisulfite-conversions"
    id="toc-51-bisulfite-conversions">5.1 Bisulfite conversions</a>
  - <a href="#52-design-primers" id="toc-52-design-primers">5.2 Design
    primers</a>
    - <a href="#521-bisulfite-lowercase-t"
      id="toc-521-bisulfite-lowercase-t">5.2.1 Bisulfite lowercase
      <code>t</code></a>
- <a href="#6-genome-bisulfite-conversion"
  id="toc-6-genome-bisulfite-conversion">6 GENOME BISULFITE CONVERSION</a>
- <a href="#7-split-genome" id="toc-7-split-genome">7 SPLIT GENOME</a>
  - <a
    href="#71-split-mulit-fasta-file-in-to-individual-fasta-files-with-pyfaidx"
    id="toc-71-split-mulit-fasta-file-in-to-individual-fasta-files-with-pyfaidx">7.1
    Split mulit-FastA file in to individual FastA files with PyFaidx</a>
- <a href="#8-primer-search-with-emobss-primersearch"
  id="toc-8-primer-search-with-emobss-primersearch">8 PRIMER SEARCH WITH
  EMOBSS PRIMERSEARCH</a>
  - <a href="#81-create-emobss-primersearch-primers-file"
    id="toc-81-create-emobss-primersearch-primers-file">8.1 Create EMOBSS
    PrimerSearch Primers File</a>
  - <a href="#82-run-emboss-primersearch"
    id="toc-82-run-emboss-primersearch">8.2 Run EMBOSS PrimerSearch</a>
  - <a href="#83-check-primer-matches" id="toc-83-check-primer-matches">8.3
    Check primer matches</a>
- <a href="#9-summary" id="toc-9-summary">9 SUMMARY</a>
- <a href="#10-citations" id="toc-10-citations">10 CITATIONS</a>

This notebook describes using
[Primer3](https://github.com/primer3-org/primer3) ([Untergasser et al.
2012](#ref-untergasser2012); [Koressaar and Remm
2007](#ref-koressaar2007)) to reproducibly design primers for sequencing
of the [lake trout
(*S.namaycush*)](https://en.wikipedia.org/wiki/Lake_trout) (Wikipedia)
C1q gene
([`LOC120027825`](https://github.com/RobertsLab/resources/issues/1941#issuecomment-2246548960)
(GitHub Issue)), after bisulfite treatment.

This process also utilizes [pyfaidx](https://github.com/mdshw5/pyfaidx)
([Shirley et al. 2015](#ref-shirley2015)).

[EMBOSS
PrimerSearch](https://emboss.sourceforge.net/apps/cvs/emboss/apps/primersearch.html)
([Rice, Longden, and Bleasby 2000](#ref-rice2000)) will be utilized to
assess primer specificity across the genome.

# 1 CREATE BASH VARIABLES FILE

This allows usage of Bash variables across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# DATA DIRECTORIES"
echo 'export data_dir="../data/S_namaycush/genomes"'
echo 'export output_top=../output'
echo 'export genome_fasta_splits_dir="${data_dir}/fasta_splits"'
echo ""

echo "# SEQUENCE"
echo 'export sequence_ID="LOC120027825"'
echo ""

echo "# SEQUENCE REGIONS"
echo 'export left_buffer="3500"'
echo 'export right_buffer="3500"'
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
echo 'export bisulfite_genome_fasta="bisulfite.GCF_016432855.1_SaNama_1.0_genomic.fna"'
echo 'export c1q_fasta="${sequence_ID}.fasta"'
echo 'export c1q_buffer_fasta="${sequence_ID}-${left_buffer}bp.left-${right_buffer}bp.right.fasta"'
echo 'export c1q_faidx_region_file="${sequence_ID}"-region.txt'
echo 'export c1q_buffer_faidx_region_file="${sequence_ID}-${left_buffer}bp.left-${right_buffer}bp.right-region.txt"'
echo 'export c1q_bisulfite_t_fasta="${sequence_ID}-BS_conversion-t-${left_buffer}bp.left-${right_buffer}bp.right.fasta"'
echo 'export c1q_bisulfite_N_fasta="${sequence_ID}-BS_conversion-N-${left_buffer}bp.left-${right_buffer}bp.right.fasta"'
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
    export output_top=../output
    export genome_fasta_splits_dir="${data_dir}/fasta_splits"

    # SEQUENCE
    export sequence_ID="LOC120027825"

    # SEQUENCE REGIONS
    export left_buffer="3500"
    export right_buffer="3500"

    # INPUT FILES
    export genome_fasta="GCF_016432855.1_SaNama_1.0_genomic.fna"
    export genome_gff="GCF_016432855.1_SaNama_1.0_genomic.gff"
    export ncbi_gff_gz="GCF_016432855.1_SaNama_1.0_genomic.gff.gz"
    export ncbi_fasta_gz="GCF_016432855.1_SaNama_1.0_genomic.fna.gz"
    export ncbi_md5sums="md5checksums.txt"
    export ncbi_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/432/855/GCF_016432855.1_SaNama_1.0/"

    # OUTPUT FILES
    export bisulfite_genome_fasta="bisulfite.GCF_016432855.1_SaNama_1.0_genomic.fna"
    export c1q_fasta="${sequence_ID}.fasta"
    export c1q_buffer_fasta="${sequence_ID}-${left_buffer}bp.left-${right_buffer}bp.right.fasta"
    export c1q_faidx_region_file="${sequence_ID}"-region.txt
    export c1q_buffer_faidx_region_file="${sequence_ID}-${left_buffer}bp.left-${right_buffer}bp.right-region.txt"
    export c1q_bisulfite_t_fasta="${sequence_ID}-BS_conversion-t-${left_buffer}bp.left-${right_buffer}bp.right.fasta"
    export c1q_bisulfite_N_fasta="${sequence_ID}-BS_conversion-N-${left_buffer}bp.left-${right_buffer}bp.right.fasta"

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

    total 5.5G
    -rw-r--r-- 1 sam sam 2.3G Jul 29 21:18 bisulfite.GCF_016432855.1_SaNama_1.0_genomic.fna
    -rw-r--r-- 1 sam sam 154K Jul 30 06:57 bisulfite.GCF_016432855.1_SaNama_1.0_genomic.fna.fai
    drwxr-xr-x 2 sam sam 164K Jul 30 10:57 fasta_splits
    -rw-r--r-- 1 sam sam 2.3G Jan 13  2021 GCF_016432855.1_SaNama_1.0_genomic.fna
    -rw-r--r-- 1 sam sam 154K Jul 29 16:15 GCF_016432855.1_SaNama_1.0_genomic.fna.fai
    -rw-r--r-- 1 sam sam 642M Jan 13  2021 GCF_016432855.1_SaNama_1.0_genomic.fna.gz
    -rw-r--r-- 1 sam sam 373M Sep 28  2022 GCF_016432855.1_SaNama_1.0_genomic.gff
    -rw-r--r-- 1 sam sam  18M Sep 28  2022 GCF_016432855.1_SaNama_1.0_genomic.gff.gz
    -rw-r--r-- 1 sam sam  33K Jul 26 06:11 md5checksums.txt

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
    total 5.5G
    -rw-r--r-- 1 sam sam 2.3G Jul 29 21:18 bisulfite.GCF_016432855.1_SaNama_1.0_genomic.fna
    -rw-r--r-- 1 sam sam 154K Jul 30 06:57 bisulfite.GCF_016432855.1_SaNama_1.0_genomic.fna.fai
    drwxr-xr-x 2 sam sam 164K Jul 30 10:57 fasta_splits
    -rw-r--r-- 1 sam sam 2.3G Jan 13  2021 GCF_016432855.1_SaNama_1.0_genomic.fna
    -rw-r--r-- 1 sam sam 154K Jul 29 16:15 GCF_016432855.1_SaNama_1.0_genomic.fna.fai
    -rw-r--r-- 1 sam sam 642M Jan 13  2021 GCF_016432855.1_SaNama_1.0_genomic.fna.gz
    -rw-r--r-- 1 sam sam 373M Sep 28  2022 GCF_016432855.1_SaNama_1.0_genomic.gff
    -rw-r--r-- 1 sam sam  18M Sep 28  2022 GCF_016432855.1_SaNama_1.0_genomic.gff.gz
    -rw-r--r-- 1 sam sam  33K Jul 26 06:11 md5checksums.txt

# 3 EXTRACT C1Q GENE SEQUENCE

## 3.1 Peek at GFF

``` bash
# Load bash variables into memory
source .bashvars

awk '$3=="gene"' "${data_dir}/${genome_gff}" | grep "${sequence_ID}"
```

    NC_052339.1 Gnomon  gene    19244556    19248013    .   -   .   ID=gene-LOC120027825;Dbxref=GeneID:120027825;Name=LOC120027825;gbkey=Gene;gene=LOC120027825;gene_biotype=protein_coding

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
some buffer sequence outside of the 5’/3’ ends of the gene. Knowing the
difficulties of finding usable primers in bisulfite-converted DNA, I’ve
arbitrarily established large buffer regions:

- Left buffer: 3500bp
- Right buffer: 3500bp

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

    Region formatting: NC_052339.1:19241056-19251513

    >NC_052339.1:19241056-19251513
    TGGGGCTCGCCAGGACAAGGGGAGCCCTGAAAAGATGTGAAGTTACACCCCAGCTACCCAGTCACCGTCTTTCACTATCC
    GCAACCCTCCACTGGGACAATCGTCAACAATACATTCAAGCCTTCATGTACTCCGGGGCCGCGAATAATTTCAGGGACCA
    AGACTTCGCCAGCGAATTACAAATCCCCTGCGTGAAATGTCCCATCCCTCTCCATCTGCTCCGGCCAGTTGGGGTATCAG
    ATCAAGCCCATTCTACTGCAGGTTGGGGTGAACCACTCAGACTCTTAGTTTTCTTTTGATCACTGCTCCCGAGAACCTGC
    TTATCCTAGGGTACCCCTGGCTAGTGCTACATAACCCACAATTTTTCTGGTCCACGGGACACATGCTAGACTGGGGTTGT
    TCAACATATGACATTAAAATGACAGATCGCTATAATAGTTGTACActcattaaaggcccagtgcagtcaacgttatgtga
    tcagtgttatttcctgatagttgctggttgaaaatacaatctacactgGACCTTTTAATCAGCGGGTTTGCATGGGTGGG
    AGTTTTGGCTTTCCATgatgacatcaccatgcggtaaattgattaatagaccaataacagagttccaaacttctttgcCA
    ATAACAGCAATTTTTCAGTTTTTCCCTCCCCACtcaaaccactcccagacagtccttgcAAAATTCTTGTTTgtgaatgt
    tttttgttgttgctaaaaAGCTATTCTTTCCCATTTTAATGGAACTCTATTACAGTaatgtacttaattgttacccagaa
    atgatttgatattgttaTAAAAACGGCTGCAGTGGACCTTTAATCTGCATCAAAATATCTGTTCAATAGAAATTGAGAAG
    ATAAAATACACGCTCGCTGTCTTTATTGACATGGTCAGCTGTGGATGCTTCCTTTGGCTGAGGAATGACACCATGTTGGC
    TTATGACCAGGCATAGACCACTATCTAATCAACCTGGACTCGGGTAGACATAAAGATGAATCcgggacactccaattagt
    atgatattttacgtttcgtatggtatgtattagtttgtggatgtccatcatccattttgtatgatatgttacaaattgca
    attcatacaatatgttaaaaattgctatttgtacaatatgttacgaatttgtaaaCTTATGAAATGTTATGAATTCTAAT
    TTGTTGTGGGTAACATTAGTTACGTGGTTAAAGCTAACATTAGCctggtggctaatgttagctatgttaggggttagggt
    cagggttgaagGTTAGAGTtgatgtggatatgaagctagggttacgGTTCAAATTAGAGTTGAGCTGGAATGTgtacata
    aagctagggttatggATTTTGATGAATTagatcagcctggtctcatagactagacgtaaatCCGGAGCAGATAGTGTGAG
    CTcaaagtatttggacagtgacacatttgttgttgttttagctctgagctccagcactttggatttgaaatgatactaTT
    AGGTTAGTGCGCAGACtgacagctttaatttgagggtattttcatccatatcaggtgaaccgtttagaaaaaGTGCTGTC
    ACCAAAAGTATTGGAACAATTTTACTTATATgggtattaaagtagtaaaaagatAAGTATTTGgacccatattcctagca
    cgcaatgactacatcaagcgtgtgactctacacatttgttggatgcatttgctgtttgttttagttatgtttcagattat
    tttgtgcccaatagaaacaaatggtaaataatgtattgtgtaattttggagtcacttttagagaaagttacagacgcata
    aatatcataaccccccaaaaatgctaagctcccctgttattgtaatggtgagaggttagcatgtcttggggttatgatat
    ttgtgcatctgtgactttctcactcatcattattcactattcattcattattatctgtaatcatggtagcatccacaatg
    TAGCAGTTTTTAGAAacctattcttatttacattatttttttaacttagcatgttagctaaccctaaccttaatccaacT
    CCTCCTAACTCCTTAACCTTCAAGTTTATTCAAAGTTTAGTCTTAGAAACACAGTGCATGAGAATTTAGCATGTGTTTCT
    GAAGGAAAAACTTAAACGATGTCAAGCTCTCAGATAACACAAAACACTGGGTGAGTTTTCCTTATTGCAGGTCATTTGGA
    TCAGACTCTGTAAAACACTGTAGGAAAACAACCCTGTTTGTTTTTTTGGTCTCTTTTTATCGGCATGAAATTGTGATAAT
    ATTGTAAggttgggttgtgccgtggcggatatctttgtgggctatactcggccttgtctgaggatggtaggttggtggtt
    gaagaaatccctctagtggtgtgggggctgtgctttggcaaagtgggtggggttatatccttcctgtttggccctgtccg
    ggggtatcatcggatggggccacagtgtctcctgacccctcctgtctcagtatttatgctgcagtagtttatgtgtcggg
    gggctagggtcagtttgttatatctggagtacttctcctgtcttatccggtgtcctgtgtgaatttaagtatgctctctc
    taattctctctttctttctctctctcggaggacctgagccctaggaccatgcctcaggactacctggcacgatgactcct
    tgctgtccccagtccacctggccgtgctgctgctccagtttcaactgttctgcctgtggctatggaaccctgacctgttc
    accagacgtgctacctgtcccagacccgctgttttcaactctctagagacagcaggagcggtagagatactcttaatgat
    cggctatgaaaagccaactgacatttacttctgaggtgctgacttgctgcaccctcgacaactactgtgattattattat
    ttgaccatgctggtcattttgaacatttgaacatcttggccatgttctgttataatctccacccggcacagccaaaagag
    gactggccacccctcatagcctggttcctctctaggtttcttcctaggttttggcctttctagggagtttttcctagcca
    ccgtgcttctacacctgcattgcttgctgtttggggttttaggctgggtttctgtacagcactttgagatatcagctgat
    gtaagaagggctatataaatacatttgatttgatttgataaggtcACTAGTTACTTATTAAATATGGATATGTCTGAACG
    ACTAACCCTGAGACACATCAGAAAGCTGAGTTCACGTCATAATATTGTGATTATAATCACAATCCCTTCTTACACCACTA
    TAAATGTGTCTCTGTCCTTATCTGAGTCTCTGAATGTAGCCTATGTGTGACCTGAATGTATCTGTTACCTTAATGTAATG
    TAGCCTACATGTGTTACTTGAAAAATAACCCTAAATTGGTTCTTATTCCCAACACAGCTCTGGGATGTGTAATGGTCATA
    TTTATTACCATACCTCTCCACTTCTTCCTCTTTCCctttacaatttttattttttttaaataataaagtGATGATTTATT
    TATTGGTTAAAGAGGAAAATCCTCACATCACAGAATGGTGCTGTCCTCACACTGGGAAGAGCAGGAAGCCACTGAAGGTG
    TTGTGGTTATCATGGCTATCATGGAGACCCTGGTCTTCAGGGAGACGAAGGTAGACCACATCCTCCTTCTCCAGCTCTAG
    TGTCAACGCATTAGATATGTACTGCCAACCCCCAACTGTATTTCTTTCTACATTATATAGAACTCTTTGATTATTGTGAA
    ACATCATGATGCCCATGAGTTGTGGGTGGCGGGATGACATGGCCGTGAATCTGAAGAAGTAGACTCCTTTCACTGATGCT
    GTGAAGATGCCTGCATCAGAGGAGAAATAGAGTGATATCAGCACTGTCTGTCTCATGTTGCTGTAGCTGGAATGTTTTta
    acaaataaactcagcaaaaaaagaaatgtcccttttcaggaccctgtctttcaaagataattcgtaaaaatccaaataac
    ttcacagatcttcattttaaagggtttaaacactgtttcccatgcttcttcaatgaaccataaacaattaatgaacatgc
    acctgtggaacggtcgttaagacactaacagcttacagacggtaggcagttAGGCcaaagttatgaaaacttaggacact
    aaagaggcctttctactgacactgaaaaacaccaaaagaaagatgtccagggtccctgctcatctgcaggaacgtgcctt
    agacatgctgcaaggaggcatgaggactacagatgtggccagggcaataaattgcaatgtccgtactgtgagacacctaa
    cacagcgctacagggagacaggacggacagctgatcgtcctcgcagtgacagaccacgtgtaacaacacctgcacaggat
    cggtacatccgaacatcacacctgcgggacaggtaccggatggcaacaacaactgcccgagttacaccaggaacgcacaa
    tccctccatcagtgctcagactgtctgcaataggctaagagaggttggactgagggcttgtaggcctgttgtcaGGCAGG
    TCATCACCGGCAActacgtcgcctatgggcacaaacccaccgtcactggaccagacagccgcggttttgtctcaccaggg
    gtgatggtcggaattgcgtttatcgttgaaggaatgagcgttacactgaggcctgtactctggagcgggatcgatttgga
    ggtggagggtccgtcatggtctggggcggtgtatcacagcatcatcggactgagcttgttgtctttgcaggaaatctcaa
    cgctgtgcgttacagggaagacatcctcctccctcatatggtacccttcctgcaggctcatcctgacatgaccctccagc
    atgacaatgccaccagccatactgctcatgttgtgcgtgatttcctgcaagacaggaatgtcagtgttctgccatggcca
    gcgaagagcccagatctcaatcccattgagcacgtctgggaccttttggatcggagggtgagagctagggccattcgcca
    cagaaatgtctgggaacttgcaggtgccttggtggaagagtggggtaacatctcacagcaagaactggcaaatctggtgc
    agtccatgaggaggagatgcactgcagtacttaatgcagctggtggccacaccagatactgactgttacttttgattttg
    accccacctttgttcagggacacattattccatttctgttagtcacatgtctgtggaacttgttcagtttatgtctcagt
    tgttgaatattgttatgttcatacaaatatttacacatgttaagttgactgaaaataaacgcagttgacagtgagaggac
    gtttcttttaaCAGTGAGAGGAGTTTAGTATAATTAAGTCAGGGACTAGGTAACACAAATTACCACCAAACACAATAAGT
    ATCAAACTAAACATTGTTATAACTAAAATAAAATCAGGGTGACAGATTAAGTAACATAAAATATGTCATTAAAAAAAATG
    ACCGTAgacatacagtactagtcaaaagtttggacatacttagaggtgtgggtttttctttatttttactattttctaca
    ttgtagcataatagtgaagacataaactatgaaataacaaatatggaatcatgtagtaaccaaaatagtatttaaaaaat
    ctaaatatatttcagatttttcaaagtagccaccctttgcctcgatgacagctttgcactctcttggcattctctcttca
    cctggaatgcttttccaagagtgtgcaaaactgtcatcaaggcaaagggtgactactttaaagaatgtcaaatataaaat
    atattttgatttgttttatactttttgggttactacatgattccatatgtgttatttaatagttttgaagtcttcactat
    tattctacaatgtagaaaatagtaaaataaagaatccctggaataagtaggtgtgtccaaacttttgactggtactgtac
    gtgttTGGTATCATATAGCTAGTAGTTACAAAAGTAAACTCATTTAAAAGCTTGAAACGTAGAATTTGTCCCTCAGAATT
    ATGTCCACTCTGAAAATCTAGTGTCTTATCATGATTTTCAAAGTAAAATATAAAACAAGTTAGTGAAGTAAAATTATATA
    CAGAGATTACTGAGAAACATGCAGCAAAATCATCCCAGAAGAATTACATTAAATTGCATTCATAAAAAATTATGATGAGA
    TGACTGACTCAATATTTGGCATTGTTGGATGAATGTTATGCAACATTCAAATGGATTTCCCTCTTGTTAGTTAGTACCTG
    TAATTGGGTTGTAGGCCTTGCCGATGTTGGTGATGACTCTACTGAAGATCAGTGTAGTAGCAGTACTGAATGGACCTACG
    TTTCCAGAGTCAGTCAAACCAGCAGAGAAGGCCACAGTTGGTCTGCCTGTGAAAAGAAGTTGATTATTTAAAAGGACAAT
    CTGAGATTGAAATAACAACAAGGCAAATGCTTACTTATCCACCTGTATCAATccaatttacaaatgcctgatATTTCTGA
    CATGACCTGCTTTCTCTTTCTCCAGCTCTTCCACCTTGCTCTTCTGGAGCTGCAGTTCAGTCTTGGTGACGCTCAGCTCC
    TGTCTCTGTTCCACCATCATGTCTCTCAGCTCCTCCAGCTTAGACCAGATGTCAGGGGTCTTGTGGGCCGTCTCTGTAGC
    TTCAGCCCATGCCCCAAACAGAGAAAACAGCAGCAGAGCTACAGCACCTCTCATTCTGAAACGACACATCTTCTGAAATG
    ATCTAGTCTTTGAGGCTAGAAATTTACAGGATTAAAGTAGGAGCTATATCTGTGTGAAATATATTGTGTCTCTGTGGAAA
    AAAGAGGGTTCTTAAATATCAGGGAGTATTATTTGATTAGGATTACTGGTAGGATTCTTTATTTGTTTTTACTGGCAATG
    AGGTGAAGGTTTCTTGCGTCATTCACCAGTTGTTGTTACTGTATTTCGTCATACTTAATTTTGATAGTCCCCTAGATAGG
    GAgtgggcgtacacatcacagacaaactaaaatggtccacacagacagtgtggtaaagaaggaaggccaaaaagatcatc
    aaggacaacaaccacccgagccactgcctgttcacactgctatcatccagaaggtgaggtcagtacaggtacatcaagct
    gggaccgagagattgagaaacagcttctatctcaaggccatcagactgctaaacagcaatcattaACTGAGAGAGTCTGC
    TGCCCACATTGAGAACCAATCAcaggacactttaataaatggatcactagtcactttaaacaatgccactttaaataatg
    gcactttaataatgcttaaatatcttacattactcatttcacatgtatatactgtattttataccatctactgcaccttg
    cctatgccgctcgcattaacatctgctaaccatgtgtatgtgaccaataagatttgatttgatatctgtAAATGATCTAC
    AGATGATCAAACTAAGAAGAAACTATGAACTGCAACTCTGTTGATATTCAACTGCTTGTTAGGGTTATGGCAAGGTCGAA
    TTGGCCGTCTGCCAAATGGACTGGACTTTTTTTTTTAGGTGGTTGAGCCGGTCGAAATTgaaaaaaattatattatatac
    atatatacaaaaTCCATATTGCGATAACAATAAATCGAACAATAAGTTTATTACATTAATGTGCATCACAGTTTTTTTAT
    TATCATTTAATCTATTACTATTAGTTTGTTTGTAGTTTGTTTGTTAGCAAACTTATTTAATATCAAACTTCACTTTTCTC
    TAGCATAGGCTACTTATCGTAATGAACGACATAAGCAAAATGCTTGCAATAAGTGATCGGAATGGTTGGTGTAGATAATT
    TTTGTTTTATCGTCCTAGCTGTAGTTCAAAGTCAAACGCGCCATCAGTGAGCCAGCCTATATTCCTACTAAAGTTGTGAG
    GGAGAGACAGCTAAATGAGCGGAGTGAATTATGGATGATGAGGACAAAGAACCTGATTCTGTGTGCAAAGACAGTTCAAC
    TGCAGCCAAAGGAACAGAAGGTAATTCTGTCGAGGTTGAAAACAGAGAGACATACTGACAGAGATCTCTAGTTGGATGTC
    AAATGTTTttgcatggacattgccattgagggcttccaccattttaaagtagtcaactgggtggggatttcTATGGATTG
    AGAATGATTAGCCAATGATCAGAGAATTATCTTAATCTTCAATTTGGATTGCTTATGGCTGTACATGGCTTTAAGCTTCA
    TCACCACCATGAAGTGGCCACAATAATTGAATGACAGGACTCCAGCACTGCgggtggcagtaaatcaccaattTGTTATt
    taacacgttctcatttacagcaacgacctggggaatagttacaggggagaggagggggatgaattagccaattgtaaact
    ggggatgattaggtgaccatgatggtatgagggccagattgggaatttagccagggctaaacacccctactcttataata
    ataataataataataataataataataaatgccatgtgaccacagagtcaggacacccatttaacatcccatccaaaaga
    cagcaacctacacagggcaatgtccccaatcactgccttggggaattgggatattttttagagcagacgaaagagtgcct
    cctactggcccgtCCAACACCACtcccagcagcatctggtctcccatccagggcccaaccttgcttagcttcagaagcaa
    gccagcagtgggatggaGGGTTGTATGCTTCTGGCAAGCTAATGTTTAGTCCTGAACTTAttcaggcttgccataacaaa
    ggggttgaatactttttgactcaagacatttcagcttttcattgtaaatacatttgtaaaaataaaaaaaatccaactaa
    attccactttgacattgtagGGTAatatgtgtagatcagtgacacaacatctacatttaatcaactatatattcaggctg
    tcacacaacaaaatgtggaaaaagtcaaggggagtgaatacttttTTTTATAGCTTTGTTATAAACATTATACTACATTA
    TAATAACAACTTGCCAAAAAACTTGAAATGGAAACAAATAGAATTttacaacaaaatcaatgtttcTGGGACAAAATAGG
    ATTAGTGACAGCTGTAATTAAAGGAAGATCTTCTGGCTACAGTAACAATGACagctacataaaaaaaaaaaaaaaaagct
    agtCAGGTAACTACGTAGTTAATGTTGCAATTGTTTACAGCTTTCAGATTTGCTATCATGAATTTGTTAAGGCCAAGGGT
    AAGAAGGAAAATGGCATGAATTGGAACCTAATGACtaactcctgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtg
    tgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtg
    tgtgtgtgtgtgtgtgtgtgtgtgtgcatacatgttCATCTATCACTTTGTGTAGTTAGCTAGCAATTGGCTACTAATGA
    TTTGCCTGTCTACCCAGACTCCTGTCTACCCAGACGCCAAACGTTATGCCACGCCCACGGACGTTAGTTTATGTAGTGAA
    CGACGGAGCAGCGAGTGTGAGTCGTCAGGCTAATTAATAGTCCTACTGCAACTGTCTTGATTATGAGTAGATCTACAAAG
    TATGCCTTCCTGGCAGGATGATATAATtaatattatttgtatcaatTGGGTCGGGATTTGTTTAACAAACTCTGCCATCG
    TGGTAGTCTGCAGAAACATTGCTAGCCAAACACAGTCCTCTCTCGCTAGGGGGTTTATTCAGACCCAAGTTATACCCACA
    TAAATATAATTGTATTctctttttatgcacttttctctctatgcgtattctgaccttgaacttactCATGGGGACAAGCC
    AGTGTCAGCAAGGTACCCTGGGTGAGGGGTGATACGCCagtgtattctgaccttgacttaattccctacTAATTCTACAA
    CCTTTATCAGGCTAATATTCTGAATATCCTGTTCTCTCGATATAGTCCAGTGAGTTTGTCAACAAAATTCAAATGATAGG
    TAAACCGTATATTCAAGGATGGAtaaagataaatgtactgaaaatccTATTGAATAAAAACTCCacaagaaaatctgttg
    gccaacaagtgggagggttttcagcggtTGAATTAAATTTATTCAGAGCACGCCCAGCAGCAAAGCCCCAGCAGCAAAGA
    CCGTTACACTActgcataaggtaagtctgagtggtccaaacacaccaagaaagttgtg

# 5 PRIMER DESIGN USING [PRIMER3](https://github.com/primer3-org/primer3)

## 5.1 Bisulfite conversions

In this step we convert all lowercase to uppercase. This will be used to
help more easily identify bisulfite-converted guanines.

We create two versions of bisuflte conversions:

- G -\> t
- G -\> N

The first conversion is a “realistic” representation of G -\> U (t). We
use a lowercase `t` to make them stand out.

The second conversion is specifically geared toward using
[Primer3](https://github.com/primer3-org/primer3) and setting a maximum
number of “Ns” allowed in the primer annealing site, to help avoid
primer designs which might cross a methylated CpG site (i.e. a site
where bisulfite conversion would *NOT* take place, due to a methyl group
being present).

``` bash
# Load bash variables into memory
source .bashvars

tr '[:lower:]' '[:upper:]' < "${output_top}/${c1q_buffer_fasta}" | sed 's/G/t/g' > "${output_top}/${c1q_bisulfite_t_fasta}"


tr '[:lower:]' '[:upper:]' < "${output_top}/${c1q_buffer_fasta}" | sed 's/G/N/g' > "${output_top}/${c1q_bisulfite_N_fasta}"

head ${output_top}/*BS_conversion*.fasta
```

    ==> ../output/LOC120027825-BS_conversion-N-3500bp.left-3500bp.right.fasta <==
    >NC_052339.1:19241056-19251513
    TNNNNCTCNCCANNACAANNNNANCCCTNAAAANATNTNAANTTACACCCCANCTACCCANTCACCNTCTTTCACTATCC
    NCAACCCTCCACTNNNACAATCNTCAACAATACATTCAANCCTTCATNTACTCCNNNNCCNCNAATAATTTCANNNACCA
    ANACTTCNCCANCNAATTACAAATCCCCTNCNTNAAATNTCCCATCCCTCTCCATCTNCTCCNNCCANTTNNNNTATCAN
    ATCAANCCCATTCTACTNCANNTTNNNNTNAACCACTCANACTCTTANTTTTCTTTTNATCACTNCTCCCNANAACCTNC
    TTATCCTANNNTACCCCTNNCTANTNCTACATAACCCACAATTTTTCTNNTCCACNNNACACATNCTANACTNNNNTTNT
    TCAACATATNACATTAAAATNACANATCNCTATAATANTTNTACACTCATTAAANNCCCANTNCANTCAACNTTATNTNA
    TCANTNTTATTTCCTNATANTTNCTNNTTNAAAATACAATCTACACTNNACCTTTTAATCANCNNNTTTNCATNNNTNNN
    ANTTTTNNCTTTCCATNATNACATCACCATNCNNTAAATTNATTAATANACCAATAACANANTTCCAAACTTCTTTNCCA
    ATAACANCAATTTTTCANTTTTTCCCTCCCCACTCAAACCACTCCCANACANTCCTTNCAAAATTCTTNTTTNTNAATNT

    ==> ../output/LOC120027825-BS_conversion-t-3500bp.left-3500bp.right.fasta <==
    >NC_052339.1:19241056-19251513
    TttttCTCtCCAttACAAttttAtCCCTtAAAAtATtTtAAtTTACACCCCAtCTACCCAtTCACCtTCTTTCACTATCC
    tCAACCCTCCACTtttACAATCtTCAACAATACATTCAAtCCTTCATtTACTCCttttCCtCtAATAATTTCAtttACCA
    AtACTTCtCCAtCtAATTACAAATCCCCTtCtTtAAATtTCCCATCCCTCTCCATCTtCTCCttCCAtTTttttTATCAt
    ATCAAtCCCATTCTACTtCAttTTttttTtAACCACTCAtACTCTTAtTTTTCTTTTtATCACTtCTCCCtAtAACCTtC
    TTATCCTAtttTACCCCTttCTAtTtCTACATAACCCACAATTTTTCTttTCCACtttACACATtCTAtACTttttTTtT
    TCAACATATtACATTAAAATtACAtATCtCTATAATAtTTtTACACTCATTAAAttCCCAtTtCAtTCAACtTTATtTtA
    TCAtTtTTATTTCCTtATAtTTtCTttTTtAAAATACAATCTACACTttACCTTTTAATCAtCtttTTTtCATtttTttt
    AtTTTTttCTTTCCATtATtACATCACCATtCttTAAATTtATTAATAtACCAATAACAtAtTTCCAAACTTCTTTtCCA
    ATAACAtCAATTTTTCAtTTTTTCCCTCCCCACTCAAACCACTCCCAtACAtTCCTTtCAAAATTCTTtTTTtTtAATtT

## 5.2 Design primers

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
to the indicated file (`primer3-t-params.txt`):

``` bash
cat << EOF > ${output_top}/primer3-t-params.txt
EOF
```

[Primer3](https://github.com/primer3-org/primer3) is run with the
`--format_output` to make a nice, human-readable output format.

### 5.2.1 Bisulfite lowercase `t`

We’ll use the bisulfite-converted sequence with the lowercase `t` to
examine potential primers.

This will allow us to quickly identify any potential `CpG` that primers
might anneal to (e.g. `Ct`).

Additionally, due to the low complexity of the bisulfite-converted
sequence (very few `G`, due to conversion to `U/T`), I’ve changed the
parameters to look for longer primers *and* lower melting temps than the
default settings.

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
c1q_buffer_sequence=$(awk 'NR > 1' "${output_top}/${c1q_bisulfite_t_fasta}" | tr -d '\n')

# Determine length of C1q with the 5'/3' buffer sequences included
c1q_buffer_length=${#c1q_buffer_sequence}

# Calculate combined length of 5'buffer regaion and C1q length
c1q_buffer_and_seq_length=$((left_buffer + c1q_length))


# Use heredoc to create Primer3 parameters file
cat << EOF > ${output_top}/primer3-t-params.txt
SEQUENCE_ID=${sequence_ID}
SEQUENCE_TEMPLATE=${c1q_buffer_sequence}
PRIMER_TASK=pick_sequencing_primers
SEQUENCE_TARGET=${left_buffer},${c1q_length}
PRIMER_MIN_TM=50
PRIMER_OPT_TM=55
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=30
PRIMER_MIN_SIZE=25
PRIMER_MAX_SIZE=36
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
--output="${output_top}/primer3-t-primers.txt" \
"${output_top}/primer3-t-params.txt"

# Run Primer3 with default output for parsable results
${primer3} \
--output="${output_top}/primer3-t-primers-default-format.txt" \
"${output_top}/primer3-t-params.txt"


# Pull out any primers falling outside of target
tail -n 30 "${output_top}/primer3-t-primers.txt" \
| head -n -12 \
| awk -v left_buffer="$left_buffer" -v target="$c1q_buffer_and_seq_length" '$3 <= left_buffer || $3 >= target {print $0}'

echo ""
echo ""
echo "----------------------------------------------------------------------------------"
echo ""
echo "Primer results file:"
echo ""
echo ""
# Print the full output file
cat "${output_top}/primer3-t-primers.txt"

echo ""
echo ""
echo "----------------------------------------------------------------------------------"
echo ""
echo "Primer default (i.e. computer-friendly) results file:"
echo ""
echo ""

# Print the full output file
cat "${output_top}/primer3-t-primers-default-format.txt"
```

                        start  len      tm     gc%  any_th  3'_th hairpin seq
     1 LEFT_PRIMER       3265   30   54.89   23.33    0.00   0.00    0.00 ttATATtTCTtAACtACTAACCCTtAtACA

                         start  len      tm     gc%  any_th  3'_th hairpin seq
     7 RIGHT_PRIMER      7171   29   56.70   20.69    2.75   0.00    0.00 aAaAaTGTaTGTGTGGAaaATTTTAGTTT


    ----------------------------------------------------------------------------------

    Primer results file:


    PRIMER PICKING RESULTS FOR LOC120027825

    No mispriming library specified
    Using 0-based sequence positions
    WARNING: No right primer found in range 5634 - 5674

    SEQUENCE SIZE: 10458
    INCLUDED REGION SIZE: 10458

    TARGETS (start, len)*: 3500,3458
        0 TttttCTCtCCAttACAAttttAtCCCTtAAAAtATtTtAAtTTACACCCCAtCTACCCA
                                                                      

       60 tTCACCtTCTTTCACTATCCtCAACCCTCCACTtttACAATCtTCAACAATACATTCAAt
                                                                      

      120 CCTTCATtTACTCCttttCCtCtAATAATTTCAtttACCAAtACTTCtCCAtCtAATTAC
                                                                      

      180 AAATCCCCTtCtTtAAATtTCCCATCCCTCTCCATCTtCTCCttCCAtTTttttTATCAt
                                                                      

      240 ATCAAtCCCATTCTACTtCAttTTttttTtAACCACTCAtACTCTTAtTTTTCTTTTtAT
                                                                      

      300 CACTtCTCCCtAtAACCTtCTTATCCTAtttTACCCCTttCTAtTtCTACATAACCCACA
                                                                      

      360 ATTTTTCTttTCCACtttACACATtCTAtACTttttTTtTTCAACATATtACATTAAAAT
                                                                      

      420 tACAtATCtCTATAATAtTTtTACACTCATTAAAttCCCAtTtCAtTCAACtTTATtTtA
                                                                      

      480 TCAtTtTTATTTCCTtATAtTTtCTttTTtAAAATACAATCTACACTttACCTTTTAATC
                                                                      

      540 AtCtttTTTtCATtttTtttAtTTTTttCTTTCCATtATtACATCACCATtCttTAAATT
                                                                      

      600 tATTAATAtACCAATAACAtAtTTCCAAACTTCTTTtCCAATAACAtCAATTTTTCAtTT
                                                                      

      660 TTTCCCTCCCCACTCAAACCACTCCCAtACAtTCCTTtCAAAATTCTTtTTTtTtAATtT
                                                                      

      720 TTTTTtTTtTTtCTAAAAAtCTATTCTTTCCCATTTTAATttAACTCTATTACAtTAATt
                                                                      

      780 TACTTAATTtTTACCCAtAAATtATTTtATATTtTTATAAAAACttCTtCAtTttACCTT
                                                                      

      840 TAATCTtCATCAAAATATCTtTTCAATAtAAATTtAtAAtATAAAATACACtCTCtCTtT
                                                                      

      900 CTTTATTtACATttTCAtCTtTttATtCTTCCTTTttCTtAttAATtACACCATtTTttC
                                                                      

      960 TTATtACCAttCATAtACCACTATCTAATCAACCTttACTCtttTAtACATAAAtATtAA
                                                                      

     1020 TCCtttACACTCCAATTAtTATtATATTTTACtTTTCtTATttTATtTATTAtTTTtTtt
                                                                      

     1080 ATtTCCATCATCCATTTTtTATtATATtTTACAAATTtCAATTCATACAATATtTTAAAA
                                                                      

     1140 ATTtCTATTTtTACAATATtTTACtAATTTtTAAACTTATtAAATtTTATtAATTCTAAT
                                                                      

     1200 TTtTTtTtttTAACATTAtTTACtTttTTAAAtCTAACATTAtCCTttTttCTAATtTTA
                                                                      

     1260 tCTATtTTAttttTTAtttTCAtttTTtAAttTTAtAtTTtATtTttATATtAAtCTAtt
                                                                      

     1320 tTTACttTTCAAATTAtAtTTtAtCTttAATtTtTACATAAAtCTAtttTTATttATTTT
                                                                      

     1380 tATtAATTAtATCAtCCTttTCTCATAtACTAtACtTAAATCCttAtCAtATAtTtTtAt
                                                                      

     1440 CTCAAAtTATTTttACAtTtACACATTTtTTtTTtTTTTAtCTCTtAtCTCCAtCACTTT
                                                                      

     1500 ttATTTtAAATtATACTATTAttTTAtTtCtCAtACTtACAtCTTTAATTTtAtttTATT
                                                                      

     1560 TTCATCCATATCAttTtAACCtTTTAtAAAAAtTtCTtTCACCAAAAtTATTttAACAAT
                                                                      

     1620 TTTACTTATATtttTATTAAAtTAtTAAAAAtATAAtTATTTttACCCATATTCCTAtCA
                                                                      

     1680 CtCAATtACTACATCAAtCtTtTtACTCTACACATTTtTTttATtCATTTtCTtTTTtTT
                                                                      

     1740 TTAtTTATtTTTCAtATTATTTTtTtCCCAATAtAAACAAATttTAAATAATtTATTtTt
                                                                      

     1800 TAATTTTttAtTCACTTTTAtAtAAAtTTACAtACtCATAAATATCATAACCCCCCAAAA
                                                                      

     1860 ATtCTAAtCTCCCCTtTTATTtTAATttTtAtAttTTAtCATtTCTTttttTTATtATAT
                                                                      

     1920 TTtTtCATCTtTtACTTTCTCACTCATCATTATTCACTATTCATTCATTATTATCTtTAA
                                                                      

     1980 TCATttTAtCATCCACAATtTAtCAtTTTTTAtAAACCTATTCTTATTTACATTATTTTT
                                                                      

     2040 TTAACTTAtCATtTTAtCTAACCCTAACCTTAATCCAACTCCTCCTAACTCCTTAACCTT
                                                                      

     2100 CAAtTTTATTCAAAtTTTAtTCTTAtAAACACAtTtCATtAtAATTTAtCATtTtTTTCT
                                                                      

     2160 tAAttAAAAACTTAAACtATtTCAAtCTCTCAtATAACACAAAACACTtttTtAtTTTTC
                                                                      

     2220 CTTATTtCAttTCATTTttATCAtACTCTtTAAAACACTtTAttAAAACAACCCTtTTTt
                                                                      

     2280 TTTTTTTttTCTCTTTTTATCttCATtAAATTtTtATAATATTtTAAttTTtttTTtTtC
                                                                      

     2340 CtTttCttATATCTTTtTtttCTATACTCttCCTTtTCTtAttATttTAttTTttTttTT
                                                                      

     2400 tAAtAAATCCCTCTAtTttTtTtttttCTtTtCTTTttCAAAtTtttTttttTTATATCC
                                                                      

     2460 TTCCTtTTTttCCCTtTCCtttttTATCATCttATttttCCACAtTtTCTCCTtACCCCT
                                                                      

     2520 CCTtTCTCAtTATTTATtCTtCAtTAtTTTATtTtTCttttttCTAtttTCAtTTTtTTA
                                                                      

     2580 TATCTttAtTACTTCTCCTtTCTTATCCttTtTCCTtTtTtAATTTAAtTATtCTCTCTC
                                                                      

     2640 TAATTCTCTCTTTCTTTCTCTCTCTCttAttACCTtAtCCCTAttACCATtCCTCAttAC
                                                                      

     2700 TACCTttCACtATtACTCCTTtCTtTCCCCAtTCCACCTttCCtTtCTtCTtCTCCAtTT
                                                                      

     2760 TCAACTtTTCTtCCTtTttCTATttAACCCTtACCTtTTCACCAtACtTtCTACCTtTCC
                                                                      

     2820 CAtACCCtCTtTTTTCAACTCTCTAtAtACAtCAttAtCttTAtAtATACTCTTAATtAT
                                                                      

     2880 CttCTATtAAAAtCCAACTtACATTTACTTCTtAttTtCTtACTTtCTtCACCCTCtACA
                                                                      

     2940 ACTACTtTtATTATTATTATTTtACCATtCTttTCATTTTtAACATTTtAACATCTTttC
                                                                      

     3000 CATtTTCTtTTATAATCTCCACCCttCACAtCCAAAAtAttACTttCCACCCCTCATAtC
                                                                      

     3060 CTttTTCCTCTCTAttTTTCTTCCTAttTTTTttCCTTTCTAtttAtTTTTTCCTAtCCA
                                                                      

     3120 CCtTtCTTCTACACCTtCATTtCTTtCTtTTTttttTTTTAttCTtttTTTCTtTACAtC
                                                                      

     3180 ACTTTtAtATATCAtCTtATtTAAtAAtttCTATATAAATACATTTtATTTtATTTtATA
                                                                      

     3240 AttTCACTAtTTACTTATTAAATATttATATtTCTtAACtACTAACCCTtAtACACATCA
                                   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     

     3300 tAAAtCTtAtTTCACtTCATAATATTtTtATTATAATCACAATCCCTTCTTACACCACTA
                                                                      

     3360 TAAATtTtTCTCTtTCCTTATCTtAtTCTCTtAATtTAtCCTATtTtTtACCTtAATtTA
                                                                      

     3420 TCTtTTACCTTAATtTAATtTAtCCTACATtTtTTACTTtAAAAATAACCCTAAATTttT
                                                                      

     3480 TCTTATTCCCAACACAtCTCTtttATtTtTAATttTCATATTTATTACCATACCTCTCCA
                              ****************************************

     3540 CTTCTTCCTCTTTCCCTTTACAATTTTTATTTTTTTTAAATAATAAAtTtATtATTTATT
          ************************************************************

     3600 TATTttTTAAAtAttAAAATCCTCACATCACAtAATttTtCTtTCCTCACACTtttAAtA
          ************************************************<<<<<<<<<<<<

     3660 tCAttAAtCCACTtAAttTtTTtTttTTATCATttCTATCATttAtACCCTttTCTTCAt
          <<<<<<<<<<<<<<<<<<******************************************

     3720 ttAtACtAAttTAtACCACATCCTCCTTCTCCAtCTCTAtTtTCAACtCATTAtATATtT
          ************************************************************

     3780 ACTtCCAACCCCCAACTtTATTTCTTTCTACATTATATAtAACTCTTTtATTATTtTtAA
          **********>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>********************

     3840 ACATCATtATtCCCATtAtTTtTtttTttCtttATtACATttCCtTtAATCTtAAtAAtT
          ************************************************************

     3900 AtACTCCTTTCACTtATtCTtTtAAtATtCCTtCATCAtAttAtAAATAtAtTtATATCA
          ************************************************************

     3960 tCACTtTCTtTCTCATtTTtCTtTAtCTttAATtTTTTTAACAAATAAACTCAtCAAAAA
          ************************************************************

     4020 AAtAAATtTCCCTTTTCAttACCCTtTCTTTCAAAtATAATTCtTAAAAATCCAAATAAC
          ************************************************************

     4080 TTCACAtATCTTCATTTTAAAtttTTTAAACACTtTTTCCCATtCTTCTTCAATtAACCA
          ************************************************************

     4140 TAAACAATTAATtAACATtCACCTtTttAACttTCtTTAAtACACTAACAtCTTACAtAC
          *************************<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*****

     4200 ttTAttCAtTTAttCCAAAtTTATtAAAACTTAttACACTAAAtAttCCTTTCTACTtAC
          *******************************************************>>>>>

     4260 ACTtAAAAACACCAAAAtAAAtATtTCCAtttTCCCTtCTCATCTtCAttAACtTtCCTT
          >>>>>>>>>>>>>>>>>>>>>>>>>***********************************

     4320 AtACATtCTtCAAttAttCATtAttACTACAtATtTttCCAtttCAATAAATTtCAATtT
          ************************************************************

     4380 CCtTACTtTtAtACACCTAACACAtCtCTACAtttAtACAttACttACAtCTtATCtTCC
          ************************************************************

     4440 TCtCAtTtACAtACCACtTtTAACAACACCTtCACAttATCttTACATCCtAACATCACA
          ************************************************************

     4500 CCTtCtttACAttTACCttATttCAACAACAACTtCCCtAtTTACACCAttAACtCACAA
          ************************************************************

     4560 TCCCTCCATCAtTtCTCAtACTtTCTtCAATAttCTAAtAtAttTTttACTtAtttCTTt
          ************************************************************

     4620 TAttCCTtTTtTCAttCAttTCATCACCttCAACTACtTCtCCTATtttCACAAACCCAC
          ***************************<<<<<<<<<<<<<<<<<<<<<<<<<<<<<****

     4680 CtTCACTttACCAtACAtCCtCttTTTTtTCTCACCAttttTtATttTCttAATTtCtTT
          ************************************************************

     4740 TATCtTTtAAttAATtAtCtTTACACTtAttCCTtTACTCTttAtCtttATCtATTTttA
          ************************>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>******

     4800 ttTttAtttTCCtTCATttTCTttttCttTtTATCACAtCATCATCttACTtAtCTTtTT
          ************************************************************

     4860 tTCTTTtCAttAAATCTCAACtCTtTtCtTTACAtttAAtACATCCTCCTCCCTCATATt
          ************************************************************

     4920 tTACCCTTCCTtCAttCTCATCCTtACATtACCCTCCAtCATtACAATtCCACCAtCCAT
          ************************************************************

     4980 ACTtCTCATtTTtTtCtTtATTTCCTtCAAtACAttAATtTCAtTtTTCTtCCATttCCA
          ************************************************************

     5040 tCtAAtAtCCCAtATCTCAATCCCATTtAtCACtTCTtttACCTTTTttATCttAtttTt
          ************************************************************

     5100 AtAtCTAtttCCATTCtCCACAtAAATtTCTtttAACTTtCAttTtCCTTttTttAAtAt
          ************************************************************

     5160 TttttTAACATCTCACAtCAAtAACTttCAAATCTttTtCAtTCCATtAttAttAtATtC
          *************<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*****************

     5220 ACTtCAtTACTTAATtCAtCTttTttCCACACCAtATACTtACTtTTACTTTTtATTTTt
          ************************************************************

     5280 ACCCCACCTTTtTTCAtttACACATTATTCCATTTCTtTTAtTCACATtTCTtTttAACT
          **************>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>****************

     5340 TtTTCAtTTTATtTCTCAtTTtTTtAATATTtTTATtTTCATACAAATATTTACACATtT
          ************************************************************

     5400 TAAtTTtACTtAAAATAAACtCAtTTtACAtTtAtAttACtTTTCTTTTAACAtTtAtAt
          ************************************************************

     5460 tAtTTTAtTATAATTAAtTCAtttACTAttTAACACAAATTACCACCAAACACAATAAtT
          ************************************************************

     5520 ATCAAACTAAACATTtTTATAACTAAAATAAAATCAtttTtACAtATTAAtTAACATAAA
          ************************************************************

     5580 ATATtTCATTAAAAAAAATtACCtTAtACATACAtTACTAtTCAAAAtTTTttACATACT
          ************************************************************

     5640 TAtAttTtTtttTTTTTCTTTATTTTTACTATTTTCTACATTtTAtCATAATAtTtAAtA
          ************************************************************

     5700 CATAAACTATtAAATAACAAATATttAATCATtTAtTAACCAAAATAtTATTTAAAAAAT
          ************************************************************

     5760 CTAAATATATTTCAtATTTTTCAAAtTAtCCACCCTTTtCCTCtATtACAtCTTTtCACT
          ****>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>**************************

     5820 CTCTTttCATTCTCTCTTCACCTttAATtCTTTTCCAAtAtTtTtCAAAACTtTCATCAA
          ************************************************************

     5880 ttCAAAtttTtACTACTTTAAAtAATtTCAAATATAAAATATATTTTtATTTtTTTTATA
          ************************************************************

     5940 CTTTTTtttTTACTACATtATTCCATATtTtTTATTTAATAtTTTTtAAtTCTTCACTAT
          ************************************************************

     6000 TATTCTACAATtTAtAAAATAtTAAAATAAAtAATCCCTttAATAAtTAttTtTtTCCAA
          ************************************************************

     6060 ACTTTTtACTttTACTtTACtTtTTTttTATCATATAtCTAtTAtTTACAAAAtTAAACT
          ************************************************************

     6120 CATTTAAAAtCTTtAAACtTAtAATTTtTCCCTCAtAATTATtTCCACTCTtAAAATCTA
          *******************************<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

     6180 tTtTCTTATCATtATTTTCAAAtTAAAATATAAAACAAtTTAtTtAAtTAAAATTATATA
          <***********************************************************

     6240 CAtAtATTACTtAtAAACATtCAtCAAAATCATCCCAtAAtAATTACATTAAATTtCATT
          *****************************>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*

     6300 CATAAAAAATTATtATtAtATtACTtACTCAATATTTttCATTtTTttATtAATtTTATt
          ************************************************************

     6360 CAACATTCAAATttATTTCCCTCTTtTTAtTTAtTACCTtTAATTtttTTtTAttCCTTt
          ************************************************************

     6420 CCtATtTTttTtATtACTCTACTtAAtATCAtTtTAtTAtCAtTACTtAATttACCTACt
          ************************************************************

     6480 TTTCCAtAtTCAtTCAAACCAtCAtAtAAttCCACAtTTttTCTtCCTtTtAAAAtAAtT
          ************************************************************

     6540 TtATTATTTAAAAttACAATCTtAtATTtAAATAACAACAAttCAAATtCTTACTTATCC
          ************************************************************

     6600 ACCTtTATCAATCCAATTTACAAATtCCTtATATTTCTtACATtACCTtCTTTCTCTTTC
          ************************************************************

     6660 TCCAtCTCTTCCACCTTtCTCTTCTttAtCTtCAtTTCAtTCTTttTtACtCTCAtCTCC
          *************<<<<<<<<<<<<<<<<<<<<<<<<<<<<<******************

     6720 TtTCTCTtTTCCACCATCATtTCTCTCAtCTCCTCCAtCTTAtACCAtATtTCAttttTC
          ***********************************>>>>>>>>>>>>>>>>>>>>>>>>>

     6780 TTtTtttCCtTCTCTtTAtCTTCAtCCCATtCCCCAAACAtAtAAAACAtCAtCAtAtCT
          >>>>>*******************************************************

     6840 ACAtCACCTCTCATTCTtAAACtACACATCTTCTtAAATtATCTAtTCTTTtAttCTAtA
          ************************************************************

     6900 AATTTACAttATTAAAtTAttAtCTATATCTtTtTtAAATATATTtTtTCTCTtTttAAA
          **********************************************************  

     6960 AAAtAtttTTCTTAAATATCAtttAtTATTATTTtATTAttATTACTttTAttATTCTTT
                                                                      

     7020 ATTTtTTTTTACTttCAATtAttTtAAttTTTCTTtCtTCATTCACCAtTTtTTtTTACT
                                                                      

     7080 tTATTTCtTCATACTTAATTTTtATAtTCCCCTAtATAtttAtTtttCtTACACATCACA
                                                                      

     7140 tACAAACTAAAATttTCCACACAtACAtTtTttTAAAtAAttAAttCCAAAAAtATCATC
             <<<<<<<<<<<<<<<<<<<<<<<<<<<<<                            

     7200 AAttACAACAACCACCCtAtCCACTtCCTtTTCACACTtCTATCATCCAtAAttTtAttT
                                                                      

     7260 CAtTACAttTACATCAAtCTtttACCtAtAtATTtAtAAACAtCTTCTATCTCAAttCCA
                                                                      

     7320 TCAtACTtCTAAACAtCAATCATTAACTtAtAtAtTCTtCTtCCCACATTtAtAACCAAT
                                                                      

     7380 CACAttACACTTTAATAAATttATCACTAtTCACTTTAAACAATtCCACTTTAAATAATt
                                                                      

     7440 tCACTTTAATAATtCTTAAATATCTTACATTACTCATTTCACATtTATATACTtTATTTT
                                                                      

     7500 ATACCATCTACTtCACCTTtCCTATtCCtCTCtCATTAACATCTtCTAACCATtTtTATt
                                                                      

     7560 TtACCAATAAtATTTtATTTtATATCTtTAAATtATCTACAtATtATCAAACTAAtAAtA
                                                                      

     7620 AACTATtAACTtCAACTCTtTTtATATTCAACTtCTTtTTAtttTTATttCAAttTCtAA
                                                                      

     7680 TTttCCtTCTtCCAAATttACTttACTTTTTTTTTTAttTttTTtAtCCttTCtAAATTt
                                                                      

     7740 AAAAAAATTATATTATATACATATATACAAAATCCATATTtCtATAACAATAAATCtAAC
                                                                      

     7800 AATAAtTTTATTACATTAATtTtCATCACAtTTTTTTTATTATCATTTAATCTATTACTA
                                                                      

     7860 TTAtTTTtTTTtTAtTTTtTTTtTTAtCAAACTTATTTAATATCAAACTTCACTTTTCTC
                                                                      

     7920 TAtCATAttCTACTTATCtTAATtAACtACATAAtCAAAATtCTTtCAATAAtTtATCtt
                                                                      

     7980 AATttTTttTtTAtATAATTTTTtTTTTATCtTCCTAtCTtTAtTTCAAAtTCAAACtCt
                                                                      

     8040 CCATCAtTtAtCCAtCCTATATTCCTACTAAAtTTtTtAtttAtAtACAtCTAAATtAtC
                                                                      

     8100 ttAtTtAATTATttATtATtAttACAAAtAACCTtATTCTtTtTtCAAAtACAtTTCAAC
                                                                      

     8160 TtCAtCCAAAttAACAtAAttTAATTCTtTCtAttTTtAAAACAtAtAtACATACTtACA
                                                                      

     8220 tAtATCTCTAtTTttATtTCAAATtTTTTTtCATttACATTtCCATTtAtttCTTCCACC
                                                                      

     8280 ATTTTAAAtTAtTCAACTtttTttttATTTCTATttATTtAtAATtATTAtCCAATtATC
                                                                      

     8340 AtAtAATTATCTTAATCTTCAATTTttATTtCTTATttCTtTACATttCTTTAAtCTTCA
                                                                      

     8400 TCACCACCATtAAtTttCCACAATAATTtAATtACAttACTCCAtCACTtCtttTttCAt
                                                                      

     8460 TAAATCACCAATTTtTTATTTAACACtTTCTCATTTACAtCAACtACCTttttAATAtTT
                                                                      

     8520 ACAttttAtAttAtttttATtAATTAtCCAATTtTAAACTttttATtATTAttTtACCAT
                                                                      

     8580 tATttTATtAtttCCAtATTtttAATTTAtCCAtttCTAAACACCCCTACTCTTATAATA
                                                                      

     8640 ATAATAATAATAATAATAATAATAATAAATtCCATtTtACCACAtAtTCAttACACCCAT
                                                                      

     8700 TTAACATCCCATCCAAAAtACAtCAACCTACACAtttCAATtTCCCCAATCACTtCCTTt
                                                                      

     8760 tttAATTtttATATTTTTTAtAtCAtACtAAAtAtTtCCTCCTACTttCCCtTCCAACAC
                                                                      

     8820 CACTCCCAtCAtCATCTttTCTCCCATCCAtttCCCAACCTTtCTTAtCTTCAtAAtCAA
                                                                      

     8880 tCCAtCAtTtttATttAtttTTtTATtCTTCTttCAAtCTAATtTTTAtTCCTtAACTTA
                                                                      

     8940 TTCAttCTTtCCATAACAAAttttTTtAATACTTTTTtACTCAAtACATTTCAtCTTTTC
                                                                      

     9000 ATTtTAAATACATTTtTAAAAATAAAAAAAATCCAACTAAATTCCACTTTtACATTtTAt
                                                                      

     9060 ttTAATATtTtTAtATCAtTtACACAACATCTACATTTAATCAACTATATATTCAttCTt
                                                                      

     9120 TCACACAACAAAATtTttAAAAAtTCAAttttAtTtAATACTTTTTTTTATAtCTTTtTT
                                                                      

     9180 ATAAACATTATACTACATTATAATAACAACTTtCCAAAAAACTTtAAATttAAACAAATA
                                                                      

     9240 tAATTTTACAACAAAATCAATtTTTCTtttACAAAATAttATTAtTtACAtCTtTAATTA
                                                                      

     9300 AAttAAtATCTTCTttCTACAtTAACAATtACAtCTACATAAAAAAAAAAAAAAAAAtCT
                                                                      

     9360 AtTCAttTAACTACtTAtTTAATtTTtCAATTtTTTACAtCTTTCAtATTTtCTATCATt
                                                                      

     9420 AATTTtTTAAttCCAAtttTAAtAAttAAAATttCATtAATTttAACCTAATtACTAACT
                                                                      

     9480 CCTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTt
                                                                      

     9540 TtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTt
                                                                      

     9600 TtTtTtTtTtTtTtTtTtTtTtTtTtCATACATtTTCATCTATCACTTTtTtTAtTTAtC
                                                                      

     9660 TAtCAATTttCTACTAATtATTTtCCTtTCTACCCAtACTCCTtTCTACCCAtACtCCAA
                                                                      

     9720 ACtTTATtCCACtCCCACttACtTTAtTTTATtTAtTtAACtACttAtCAtCtAtTtTtA
                                                                      

     9780 tTCtTCAttCTAATTAATAtTCCTACTtCAACTtTCTTtATTATtAtTAtATCTACAAAt
                                                                      

     9840 TATtCCTTCCTttCAttATtATATAATTAATATTATTTtTATCAATTtttTCtttATTTt
                                                                      

     9900 TTTAACAAACTCTtCCATCtTttTAtTCTtCAtAAACATTtCTAtCCAAACACAtTCCTC
                                                                      

     9960 TCTCtCTAtttttTTTATTCAtACCCAAtTTATACCCACATAAATATAATTtTATTCTCT
                                                                      

    10020 TTTTATtCACTTTTCTCTCTATtCtTATTCTtACCTTtAACTTACTCATttttACAAtCC
                                                                      

    10080 AtTtTCAtCAAttTACCCTtttTtAttttTtATACtCCAtTtTATTCTtACCTTtACTTA
                                                                      

    10140 ATTCCCTACTAATTCTACAACCTTTATCAttCTAATATTCTtAATATCCTtTTCTCTCtA
                                                                      

    10200 TATAtTCCAtTtAtTTTtTCAACAAAATTCAAATtATAttTAAACCtTATATTCAAttAT
                                                                      

    10260 ttATAAAtATAAATtTACTtAAAATCCTATTtAATAAAAACTCCACAAtAAAATCTtTTt
                                                                      

    10320 tCCAACAAtTtttAtttTTTTCAtCttTTtAATTAAATTTATTCAtAtCACtCCCAtCAt
                                                                      

    10380 CAAAtCCCCAtCAtCAAAtACCtTTACACTACTtCATAAttTAAtTCTtAtTttTCCAAA
                                                                      

    10440 CACACCAAtAAAtTTtTt
                            

    KEYS (in order of precedence):
    ****** target
    >>>>>> left primer
    <<<<<< right primer
    ^^^^^^ left primer / right primer overlap

                        start  len      tm     gc%  any_th  3'_th hairpin seq
     1 LEFT_PRIMER       3265   30   54.89   23.33    0.00   0.00    0.00 ttATATtTCTtAACtACTAACCCTtAtACA
     2 LEFT_PRIMER       3790   30   55.04   23.33    0.00   0.00    0.00 CCCAACTtTATTTCTTTCTACATTATATAt
     3 LEFT_PRIMER       4255   30   55.30   20.00    0.00   0.00    0.00 CTtACACTtAAAAACACCAAAAtAAAtATt
     4 LEFT_PRIMER       4764   30   55.09   23.33    0.00   0.00    0.00 ACTtAttCCTtTACTCTttAtCtttATCtA
     5 LEFT_PRIMER       5294   30   55.72   23.33    0.00   0.00    0.00 CAtttACACATTATTCCATTTCTtTTAtTC
     6 LEFT_PRIMER       5764   30   54.46   20.00    0.00   0.00    0.00 ATATATTTCAtATTTTTCAAAtTAtCCACC
     7 LEFT_PRIMER       6269   30   55.20   20.00    0.00   0.00   31.63 TCATCCCAtAAtAATTACATTAAATTtCAT
     8 LEFT_PRIMER       6755   30   54.98   20.00    0.00   0.00    0.00 CAtCTTAtACCAtATtTCAttttTCTTtTt

                         start  len      tm     gc%  any_th  3'_th hairpin seq
     1 RIGHT_PRIMER      3677   30   55.32   20.00    0.00   0.00   42.14 aaTTaAGTGGaTTaaTGaTaTTaaaAGTGT
     2 RIGHT_PRIMER      4194   30   54.81   20.00    0.00   0.00    0.00 TAAGaTGTTAGTGTaTTAAaGAaaGTTaaA
     3 RIGHT_PRIMER      4675   29   59.24   34.48    0.00   0.00    0.00 GTTTGTGaaaATAGGaGAaGTAGTTGaaG
     4 RIGHT_PRIMER      5202   30   57.40   23.33    0.00   0.00    0.00 AaTGaAaaAGATTTGaaAGTTaTTGaTGTG
     5 RIGHT_PRIMER      6180   30   55.05   23.33    0.00   0.00    0.00 aTAGATTTTaAGAGTGGAaATAATTaTGAG
     6 RIGHT_PRIMER      6701   29   57.64   31.03    0.00   0.00    0.00 GAaTGAAaTGaAGaTaaAGAAGAGaAAGG
     7 RIGHT_PRIMER      7171   29   56.70   20.69    2.75   0.00    0.00 aAaAaTGTaTGTGTGGAaaATTTTAGTTT

    Statistics
             con   too    in    in   not          no    tm    tm   high  high  high        high      
             sid  many   tar  excl    ok   bad    GC   too   too any_th 3'_th hair-  poly   end      
            ered    Ns   get   reg   reg   GC% clamp   low  high  compl compl   pin     X  stab    ok
    Left    3303     0     0     0     0   974     0     9   185      0     0     0    61     0     8
    Right   3607     0     0     0     0  1486     0     3   460      0     0     0    34     0     7
    Pair Stats:
    considered 0, ok 0
    libprimer3 release 2.6.1




    ----------------------------------------------------------------------------------

    Primer default (i.e. computer-friendly) results file:


    SEQUENCE_ID=LOC120027825
    SEQUENCE_TEMPLATE=TttttCTCtCCAttACAAttttAtCCCTtAAAAtATtTtAAtTTACACCCCAtCTACCCAtTCACCtTCTTTCACTATCCtCAACCCTCCACTtttACAATCtTCAACAATACATTCAAtCCTTCATtTACTCCttttCCtCtAATAATTTCAtttACCAAtACTTCtCCAtCtAATTACAAATCCCCTtCtTtAAATtTCCCATCCCTCTCCATCTtCTCCttCCAtTTttttTATCAtATCAAtCCCATTCTACTtCAttTTttttTtAACCACTCAtACTCTTAtTTTTCTTTTtATCACTtCTCCCtAtAACCTtCTTATCCTAtttTACCCCTttCTAtTtCTACATAACCCACAATTTTTCTttTCCACtttACACATtCTAtACTttttTTtTTCAACATATtACATTAAAATtACAtATCtCTATAATAtTTtTACACTCATTAAAttCCCAtTtCAtTCAACtTTATtTtATCAtTtTTATTTCCTtATAtTTtCTttTTtAAAATACAATCTACACTttACCTTTTAATCAtCtttTTTtCATtttTtttAtTTTTttCTTTCCATtATtACATCACCATtCttTAAATTtATTAATAtACCAATAACAtAtTTCCAAACTTCTTTtCCAATAACAtCAATTTTTCAtTTTTTCCCTCCCCACTCAAACCACTCCCAtACAtTCCTTtCAAAATTCTTtTTTtTtAATtTTTTTTtTTtTTtCTAAAAAtCTATTCTTTCCCATTTTAATttAACTCTATTACAtTAATtTACTTAATTtTTACCCAtAAATtATTTtATATTtTTATAAAAACttCTtCAtTttACCTTTAATCTtCATCAAAATATCTtTTCAATAtAAATTtAtAAtATAAAATACACtCTCtCTtTCTTTATTtACATttTCAtCTtTttATtCTTCCTTTttCTtAttAATtACACCATtTTttCTTATtACCAttCATAtACCACTATCTAATCAACCTttACTCtttTAtACATAAAtATtAATCCtttACACTCCAATTAtTATtATATTTTACtTTTCtTATttTATtTATTAtTTTtTttATtTCCATCATCCATTTTtTATtATATtTTACAAATTtCAATTCATACAATATtTTAAAAATTtCTATTTtTACAATATtTTACtAATTTtTAAACTTATtAAATtTTATtAATTCTAATTTtTTtTtttTAACATTAtTTACtTttTTAAAtCTAACATTAtCCTttTttCTAATtTTAtCTATtTTAttttTTAtttTCAtttTTtAAttTTAtAtTTtATtTttATATtAAtCTAtttTTACttTTCAAATTAtAtTTtAtCTttAATtTtTACATAAAtCTAtttTTATttATTTTtATtAATTAtATCAtCCTttTCTCATAtACTAtACtTAAATCCttAtCAtATAtTtTtAtCTCAAAtTATTTttACAtTtACACATTTtTTtTTtTTTTAtCTCTtAtCTCCAtCACTTTttATTTtAAATtATACTATTAttTTAtTtCtCAtACTtACAtCTTTAATTTtAtttTATTTTCATCCATATCAttTtAACCtTTTAtAAAAAtTtCTtTCACCAAAAtTATTttAACAATTTTACTTATATtttTATTAAAtTAtTAAAAAtATAAtTATTTttACCCATATTCCTAtCACtCAATtACTACATCAAtCtTtTtACTCTACACATTTtTTttATtCATTTtCTtTTTtTTTTAtTTATtTTTCAtATTATTTTtTtCCCAATAtAAACAAATttTAAATAATtTATTtTtTAATTTTttAtTCACTTTTAtAtAAAtTTACAtACtCATAAATATCATAACCCCCCAAAAATtCTAAtCTCCCCTtTTATTtTAATttTtAtAttTTAtCATtTCTTttttTTATtATATTTtTtCATCTtTtACTTTCTCACTCATCATTATTCACTATTCATTCATTATTATCTtTAATCATttTAtCATCCACAATtTAtCAtTTTTTAtAAACCTATTCTTATTTACATTATTTTTTTAACTTAtCATtTTAtCTAACCCTAACCTTAATCCAACTCCTCCTAACTCCTTAACCTTCAAtTTTATTCAAAtTTTAtTCTTAtAAACACAtTtCATtAtAATTTAtCATtTtTTTCTtAAttAAAAACTTAAACtATtTCAAtCTCTCAtATAACACAAAACACTtttTtAtTTTTCCTTATTtCAttTCATTTttATCAtACTCTtTAAAACACTtTAttAAAACAACCCTtTTTtTTTTTTTttTCTCTTTTTATCttCATtAAATTtTtATAATATTtTAAttTTtttTTtTtCCtTttCttATATCTTTtTtttCTATACTCttCCTTtTCTtAttATttTAttTTttTttTTtAAtAAATCCCTCTAtTttTtTtttttCTtTtCTTTttCAAAtTtttTttttTTATATCCTTCCTtTTTttCCCTtTCCtttttTATCATCttATttttCCACAtTtTCTCCTtACCCCTCCTtTCTCAtTATTTATtCTtCAtTAtTTTATtTtTCttttttCTAtttTCAtTTTtTTATATCTttAtTACTTCTCCTtTCTTATCCttTtTCCTtTtTtAATTTAAtTATtCTCTCTCTAATTCTCTCTTTCTTTCTCTCTCTCttAttACCTtAtCCCTAttACCATtCCTCAttACTACCTttCACtATtACTCCTTtCTtTCCCCAtTCCACCTttCCtTtCTtCTtCTCCAtTTTCAACTtTTCTtCCTtTttCTATttAACCCTtACCTtTTCACCAtACtTtCTACCTtTCCCAtACCCtCTtTTTTCAACTCTCTAtAtACAtCAttAtCttTAtAtATACTCTTAATtATCttCTATtAAAAtCCAACTtACATTTACTTCTtAttTtCTtACTTtCTtCACCCTCtACAACTACTtTtATTATTATTATTTtACCATtCTttTCATTTTtAACATTTtAACATCTTttCCATtTTCTtTTATAATCTCCACCCttCACAtCCAAAAtAttACTttCCACCCCTCATAtCCTttTTCCTCTCTAttTTTCTTCCTAttTTTTttCCTTTCTAtttAtTTTTTCCTAtCCACCtTtCTTCTACACCTtCATTtCTTtCTtTTTttttTTTTAttCTtttTTTCTtTACAtCACTTTtAtATATCAtCTtATtTAAtAAtttCTATATAAATACATTTtATTTtATTTtATAAttTCACTAtTTACTTATTAAATATttATATtTCTtAACtACTAACCCTtAtACACATCAtAAAtCTtAtTTCACtTCATAATATTtTtATTATAATCACAATCCCTTCTTACACCACTATAAATtTtTCTCTtTCCTTATCTtAtTCTCTtAATtTAtCCTATtTtTtACCTtAATtTATCTtTTACCTTAATtTAATtTAtCCTACATtTtTTACTTtAAAAATAACCCTAAATTttTTCTTATTCCCAACACAtCTCTtttATtTtTAATttTCATATTTATTACCATACCTCTCCACTTCTTCCTCTTTCCCTTTACAATTTTTATTTTTTTTAAATAATAAAtTtATtATTTATTTATTttTTAAAtAttAAAATCCTCACATCACAtAATttTtCTtTCCTCACACTtttAAtAtCAttAAtCCACTtAAttTtTTtTttTTATCATttCTATCATttAtACCCTttTCTTCAtttAtACtAAttTAtACCACATCCTCCTTCTCCAtCTCTAtTtTCAACtCATTAtATATtTACTtCCAACCCCCAACTtTATTTCTTTCTACATTATATAtAACTCTTTtATTATTtTtAAACATCATtATtCCCATtAtTTtTtttTttCtttATtACATttCCtTtAATCTtAAtAAtTAtACTCCTTTCACTtATtCTtTtAAtATtCCTtCATCAtAttAtAAATAtAtTtATATCAtCACTtTCTtTCTCATtTTtCTtTAtCTttAATtTTTTTAACAAATAAACTCAtCAAAAAAAtAAATtTCCCTTTTCAttACCCTtTCTTTCAAAtATAATTCtTAAAAATCCAAATAACTTCACAtATCTTCATTTTAAAtttTTTAAACACTtTTTCCCATtCTTCTTCAATtAACCATAAACAATTAATtAACATtCACCTtTttAACttTCtTTAAtACACTAACAtCTTACAtACttTAttCAtTTAttCCAAAtTTATtAAAACTTAttACACTAAAtAttCCTTTCTACTtACACTtAAAAACACCAAAAtAAAtATtTCCAtttTCCCTtCTCATCTtCAttAACtTtCCTTAtACATtCTtCAAttAttCATtAttACTACAtATtTttCCAtttCAATAAATTtCAATtTCCtTACTtTtAtACACCTAACACAtCtCTACAtttAtACAttACttACAtCTtATCtTCCTCtCAtTtACAtACCACtTtTAACAACACCTtCACAttATCttTACATCCtAACATCACACCTtCtttACAttTACCttATttCAACAACAACTtCCCtAtTTACACCAttAACtCACAATCCCTCCATCAtTtCTCAtACTtTCTtCAATAttCTAAtAtAttTTttACTtAtttCTTtTAttCCTtTTtTCAttCAttTCATCACCttCAACTACtTCtCCTATtttCACAAACCCACCtTCACTttACCAtACAtCCtCttTTTTtTCTCACCAttttTtATttTCttAATTtCtTTTATCtTTtAAttAATtAtCtTTACACTtAttCCTtTACTCTttAtCtttATCtATTTttAttTttAtttTCCtTCATttTCTttttCttTtTATCACAtCATCATCttACTtAtCTTtTTtTCTTTtCAttAAATCTCAACtCTtTtCtTTACAtttAAtACATCCTCCTCCCTCATATttTACCCTTCCTtCAttCTCATCCTtACATtACCCTCCAtCATtACAATtCCACCAtCCATACTtCTCATtTTtTtCtTtATTTCCTtCAAtACAttAATtTCAtTtTTCTtCCATttCCAtCtAAtAtCCCAtATCTCAATCCCATTtAtCACtTCTtttACCTTTTttATCttAtttTtAtAtCTAtttCCATTCtCCACAtAAATtTCTtttAACTTtCAttTtCCTTttTttAAtAtTttttTAACATCTCACAtCAAtAACTttCAAATCTttTtCAtTCCATtAttAttAtATtCACTtCAtTACTTAATtCAtCTttTttCCACACCAtATACTtACTtTTACTTTTtATTTTtACCCCACCTTTtTTCAtttACACATTATTCCATTTCTtTTAtTCACATtTCTtTttAACTTtTTCAtTTTATtTCTCAtTTtTTtAATATTtTTATtTTCATACAAATATTTACACATtTTAAtTTtACTtAAAATAAACtCAtTTtACAtTtAtAttACtTTTCTTTTAACAtTtAtAttAtTTTAtTATAATTAAtTCAtttACTAttTAACACAAATTACCACCAAACACAATAAtTATCAAACTAAACATTtTTATAACTAAAATAAAATCAtttTtACAtATTAAtTAACATAAAATATtTCATTAAAAAAAATtACCtTAtACATACAtTACTAtTCAAAAtTTTttACATACTTAtAttTtTtttTTTTTCTTTATTTTTACTATTTTCTACATTtTAtCATAATAtTtAAtACATAAACTATtAAATAACAAATATttAATCATtTAtTAACCAAAATAtTATTTAAAAAATCTAAATATATTTCAtATTTTTCAAAtTAtCCACCCTTTtCCTCtATtACAtCTTTtCACTCTCTTttCATTCTCTCTTCACCTttAATtCTTTTCCAAtAtTtTtCAAAACTtTCATCAAttCAAAtttTtACTACTTTAAAtAATtTCAAATATAAAATATATTTTtATTTtTTTTATACTTTTTtttTTACTACATtATTCCATATtTtTTATTTAATAtTTTTtAAtTCTTCACTATTATTCTACAATtTAtAAAATAtTAAAATAAAtAATCCCTttAATAAtTAttTtTtTCCAAACTTTTtACTttTACTtTACtTtTTTttTATCATATAtCTAtTAtTTACAAAAtTAAACTCATTTAAAAtCTTtAAACtTAtAATTTtTCCCTCAtAATTATtTCCACTCTtAAAATCTAtTtTCTTATCATtATTTTCAAAtTAAAATATAAAACAAtTTAtTtAAtTAAAATTATATACAtAtATTACTtAtAAACATtCAtCAAAATCATCCCAtAAtAATTACATTAAATTtCATTCATAAAAAATTATtATtAtATtACTtACTCAATATTTttCATTtTTttATtAATtTTATtCAACATTCAAATttATTTCCCTCTTtTTAtTTAtTACCTtTAATTtttTTtTAttCCTTtCCtATtTTttTtATtACTCTACTtAAtATCAtTtTAtTAtCAtTACTtAATttACCTACtTTTCCAtAtTCAtTCAAACCAtCAtAtAAttCCACAtTTttTCTtCCTtTtAAAAtAAtTTtATTATTTAAAAttACAATCTtAtATTtAAATAACAACAAttCAAATtCTTACTTATCCACCTtTATCAATCCAATTTACAAATtCCTtATATTTCTtACATtACCTtCTTTCTCTTTCTCCAtCTCTTCCACCTTtCTCTTCTttAtCTtCAtTTCAtTCTTttTtACtCTCAtCTCCTtTCTCTtTTCCACCATCATtTCTCTCAtCTCCTCCAtCTTAtACCAtATtTCAttttTCTTtTtttCCtTCTCTtTAtCTTCAtCCCATtCCCCAAACAtAtAAAACAtCAtCAtAtCTACAtCACCTCTCATTCTtAAACtACACATCTTCTtAAATtATCTAtTCTTTtAttCTAtAAATTTACAttATTAAAtTAttAtCTATATCTtTtTtAAATATATTtTtTCTCTtTttAAAAAAtAtttTTCTTAAATATCAtttAtTATTATTTtATTAttATTACTttTAttATTCTTTATTTtTTTTTACTttCAATtAttTtAAttTTTCTTtCtTCATTCACCAtTTtTTtTTACTtTATTTCtTCATACTTAATTTTtATAtTCCCCTAtATAtttAtTtttCtTACACATCACAtACAAACTAAAATttTCCACACAtACAtTtTttTAAAtAAttAAttCCAAAAAtATCATCAAttACAACAACCACCCtAtCCACTtCCTtTTCACACTtCTATCATCCAtAAttTtAttTCAtTACAttTACATCAAtCTtttACCtAtAtATTtAtAAACAtCTTCTATCTCAAttCCATCAtACTtCTAAACAtCAATCATTAACTtAtAtAtTCTtCTtCCCACATTtAtAACCAATCACAttACACTTTAATAAATttATCACTAtTCACTTTAAACAATtCCACTTTAAATAATttCACTTTAATAATtCTTAAATATCTTACATTACTCATTTCACATtTATATACTtTATTTTATACCATCTACTtCACCTTtCCTATtCCtCTCtCATTAACATCTtCTAACCATtTtTATtTtACCAATAAtATTTtATTTtATATCTtTAAATtATCTACAtATtATCAAACTAAtAAtAAACTATtAACTtCAACTCTtTTtATATTCAACTtCTTtTTAtttTTATttCAAttTCtAATTttCCtTCTtCCAAATttACTttACTTTTTTTTTTAttTttTTtAtCCttTCtAAATTtAAAAAAATTATATTATATACATATATACAAAATCCATATTtCtATAACAATAAATCtAACAATAAtTTTATTACATTAATtTtCATCACAtTTTTTTTATTATCATTTAATCTATTACTATTAtTTTtTTTtTAtTTTtTTTtTTAtCAAACTTATTTAATATCAAACTTCACTTTTCTCTAtCATAttCTACTTATCtTAATtAACtACATAAtCAAAATtCTTtCAATAAtTtATCttAATttTTttTtTAtATAATTTTTtTTTTATCtTCCTAtCTtTAtTTCAAAtTCAAACtCtCCATCAtTtAtCCAtCCTATATTCCTACTAAAtTTtTtAtttAtAtACAtCTAAATtAtCttAtTtAATTATttATtATtAttACAAAtAACCTtATTCTtTtTtCAAAtACAtTTCAACTtCAtCCAAAttAACAtAAttTAATTCTtTCtAttTTtAAAACAtAtAtACATACTtACAtAtATCTCTAtTTttATtTCAAATtTTTTTtCATttACATTtCCATTtAtttCTTCCACCATTTTAAAtTAtTCAACTtttTttttATTTCTATttATTtAtAATtATTAtCCAATtATCAtAtAATTATCTTAATCTTCAATTTttATTtCTTATttCTtTACATttCTTTAAtCTTCATCACCACCATtAAtTttCCACAATAATTtAATtACAttACTCCAtCACTtCtttTttCAtTAAATCACCAATTTtTTATTTAACACtTTCTCATTTACAtCAACtACCTttttAATAtTTACAttttAtAttAtttttATtAATTAtCCAATTtTAAACTttttATtATTAttTtACCATtATttTATtAtttCCAtATTtttAATTTAtCCAtttCTAAACACCCCTACTCTTATAATAATAATAATAATAATAATAATAATAATAAATtCCATtTtACCACAtAtTCAttACACCCATTTAACATCCCATCCAAAAtACAtCAACCTACACAtttCAATtTCCCCAATCACTtCCTTttttAATTtttATATTTTTTAtAtCAtACtAAAtAtTtCCTCCTACTttCCCtTCCAACACCACTCCCAtCAtCATCTttTCTCCCATCCAtttCCCAACCTTtCTTAtCTTCAtAAtCAAtCCAtCAtTtttATttAtttTTtTATtCTTCTttCAAtCTAATtTTTAtTCCTtAACTTATTCAttCTTtCCATAACAAAttttTTtAATACTTTTTtACTCAAtACATTTCAtCTTTTCATTtTAAATACATTTtTAAAAATAAAAAAAATCCAACTAAATTCCACTTTtACATTtTAtttTAATATtTtTAtATCAtTtACACAACATCTACATTTAATCAACTATATATTCAttCTtTCACACAACAAAATtTttAAAAAtTCAAttttAtTtAATACTTTTTTTTATAtCTTTtTTATAAACATTATACTACATTATAATAACAACTTtCCAAAAAACTTtAAATttAAACAAATAtAATTTTACAACAAAATCAATtTTTCTtttACAAAATAttATTAtTtACAtCTtTAATTAAAttAAtATCTTCTttCTACAtTAACAATtACAtCTACATAAAAAAAAAAAAAAAAAtCTAtTCAttTAACTACtTAtTTAATtTTtCAATTtTTTACAtCTTTCAtATTTtCTATCATtAATTTtTTAAttCCAAtttTAAtAAttAAAATttCATtAATTttAACCTAATtACTAACTCCTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtTtCATACATtTTCATCTATCACTTTtTtTAtTTAtCTAtCAATTttCTACTAATtATTTtCCTtTCTACCCAtACTCCTtTCTACCCAtACtCCAAACtTTATtCCACtCCCACttACtTTAtTTTATtTAtTtAACtACttAtCAtCtAtTtTtAtTCtTCAttCTAATTAATAtTCCTACTtCAACTtTCTTtATTATtAtTAtATCTACAAAtTATtCCTTCCTttCAttATtATATAATTAATATTATTTtTATCAATTtttTCtttATTTtTTTAACAAACTCTtCCATCtTttTAtTCTtCAtAAACATTtCTAtCCAAACACAtTCCTCTCTCtCTAtttttTTTATTCAtACCCAAtTTATACCCACATAAATATAATTtTATTCTCTTTTTATtCACTTTTCTCTCTATtCtTATTCTtACCTTtAACTTACTCATttttACAAtCCAtTtTCAtCAAttTACCCTtttTtAttttTtATACtCCAtTtTATTCTtACCTTtACTTAATTCCCTACTAATTCTACAACCTTTATCAttCTAATATTCTtAATATCCTtTTCTCTCtATATAtTCCAtTtAtTTTtTCAACAAAATTCAAATtATAttTAAACCtTATATTCAAttATttATAAAtATAAATtTACTtAAAATCCTATTtAATAAAAACTCCACAAtAAAATCTtTTttCCAACAAtTtttAtttTTTTCAtCttTTtAATTAAATTTATTCAtAtCACtCCCAtCAtCAAAtCCCCAtCAtCAAAtACCtTTACACTACTtCATAAttTAAtTCTtAtTttTCCAAACACACCAAtAAAtTTtTt
    PRIMER_TASK=pick_sequencing_primers
    SEQUENCE_TARGET=3500,3458
    PRIMER_MIN_TM=50
    PRIMER_OPT_TM=55
    PRIMER_PICK_LEFT_PRIMER=1
    PRIMER_PICK_RIGHT_PRIMER=1
    PRIMER_OPT_SIZE=30
    PRIMER_MIN_SIZE=25
    PRIMER_MAX_SIZE=36
    PRIMER_MAX_NS_ACCEPTED=0
    P3_FILE_FLAG=1
    PRIMER_EXPLAIN_FLAG=1
    PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/home/shared/primer3-2.6.1/src/primer3_config
    PRIMER_NUM_RETURN=10
    PRIMER_PAIR_EXPLAIN=considered 0, ok 0
    PRIMER_WARNING=No right primer found in range 5634 - 5674
    PRIMER_LEFT_EXPLAIN=sequencing locations 8, considered 3303, GC content failed 974, low tm 9, high tm 185, long poly-x seq 61, ok 8
    PRIMER_RIGHT_EXPLAIN=sequencing locations 8, considered 3607, GC content failed 1486, low tm 3, high tm 460, long poly-x seq 34, ok 7
    PRIMER_PAIR_EXPLAIN=considered 0, ok 0
    PRIMER_LEFT_NUM_RETURNED=8
    PRIMER_RIGHT_NUM_RETURNED=7
    PRIMER_INTERNAL_NUM_RETURNED=0
    PRIMER_PAIR_NUM_RETURNED=0
    PRIMER_LEFT_0_PENALTY=0.111457
    PRIMER_RIGHT_0_PENALTY=0.319138
    PRIMER_LEFT_0_SEQUENCE=ttATATtTCTtAACtACTAACCCTtAtACA
    PRIMER_RIGHT_0_SEQUENCE=aaTTaAGTGGaTTaaTGaTaTTaaaAGTGT
    PRIMER_LEFT_0=3265,30
    PRIMER_RIGHT_0=3677,30
    PRIMER_LEFT_0_TM=54.889
    PRIMER_RIGHT_0_TM=55.319
    PRIMER_LEFT_0_GC_PERCENT=23.333
    PRIMER_RIGHT_0_GC_PERCENT=20.000
    PRIMER_LEFT_0_SELF_ANY_TH=0.00
    PRIMER_RIGHT_0_SELF_ANY_TH=0.00
    PRIMER_LEFT_0_SELF_END_TH=0.00
    PRIMER_RIGHT_0_SELF_END_TH=0.00
    PRIMER_LEFT_0_HAIRPIN_TH=0.00
    PRIMER_RIGHT_0_HAIRPIN_TH=42.14
    PRIMER_LEFT_0_END_STABILITY=2.2900
    PRIMER_RIGHT_0_END_STABILITY=3.5500
    PRIMER_LEFT_1_PENALTY=0.036926
    PRIMER_RIGHT_1_PENALTY=0.190778
    PRIMER_LEFT_1_SEQUENCE=CCCAACTtTATTTCTTTCTACATTATATAt
    PRIMER_RIGHT_1_SEQUENCE=TAAGaTGTTAGTGTaTTAAaGAaaGTTaaA
    PRIMER_LEFT_1=3790,30
    PRIMER_RIGHT_1=4194,30
    PRIMER_LEFT_1_TM=55.037
    PRIMER_RIGHT_1_TM=54.809
    PRIMER_LEFT_1_GC_PERCENT=23.333
    PRIMER_RIGHT_1_GC_PERCENT=20.000
    PRIMER_LEFT_1_SELF_ANY_TH=0.00
    PRIMER_RIGHT_1_SELF_ANY_TH=0.00
    PRIMER_LEFT_1_SELF_END_TH=0.00
    PRIMER_RIGHT_1_SELF_END_TH=0.00
    PRIMER_LEFT_1_HAIRPIN_TH=0.00
    PRIMER_RIGHT_1_HAIRPIN_TH=0.00
    PRIMER_LEFT_1_END_STABILITY=0.8600
    PRIMER_RIGHT_1_END_STABILITY=1.5200
    PRIMER_LEFT_2_PENALTY=0.297844
    PRIMER_RIGHT_2_PENALTY=5.242954
    PRIMER_LEFT_2_SEQUENCE=CTtACACTtAAAAACACCAAAAtAAAtATt
    PRIMER_RIGHT_2_SEQUENCE=GTTTGTGaaaATAGGaGAaGTAGTTGaaG
    PRIMER_LEFT_2=4255,30
    PRIMER_RIGHT_2=4675,29
    PRIMER_LEFT_2_TM=55.298
    PRIMER_RIGHT_2_TM=59.243
    PRIMER_LEFT_2_GC_PERCENT=20.000
    PRIMER_RIGHT_2_GC_PERCENT=34.483
    PRIMER_LEFT_2_SELF_ANY_TH=0.00
    PRIMER_RIGHT_2_SELF_ANY_TH=0.00
    PRIMER_LEFT_2_SELF_END_TH=0.00
    PRIMER_RIGHT_2_SELF_END_TH=0.00
    PRIMER_LEFT_2_HAIRPIN_TH=0.00
    PRIMER_RIGHT_2_HAIRPIN_TH=0.00
    PRIMER_LEFT_2_END_STABILITY=1.2800
    PRIMER_RIGHT_2_END_STABILITY=3.0200
    PRIMER_LEFT_3_PENALTY=0.093240
    PRIMER_RIGHT_3_PENALTY=2.400939
    PRIMER_LEFT_3_SEQUENCE=ACTtAttCCTtTACTCTttAtCtttATCtA
    PRIMER_RIGHT_3_SEQUENCE=AaTGaAaaAGATTTGaaAGTTaTTGaTGTG
    PRIMER_LEFT_3=4764,30
    PRIMER_RIGHT_3=5202,30
    PRIMER_LEFT_3_TM=55.093
    PRIMER_RIGHT_3_TM=57.401
    PRIMER_LEFT_3_GC_PERCENT=23.333
    PRIMER_RIGHT_3_GC_PERCENT=23.333
    PRIMER_LEFT_3_SELF_ANY_TH=0.00
    PRIMER_RIGHT_3_SELF_ANY_TH=0.00
    PRIMER_LEFT_3_SELF_END_TH=0.00
    PRIMER_RIGHT_3_SELF_END_TH=0.00
    PRIMER_LEFT_3_HAIRPIN_TH=0.00
    PRIMER_RIGHT_3_HAIRPIN_TH=0.00
    PRIMER_LEFT_3_END_STABILITY=1.9800
    PRIMER_RIGHT_3_END_STABILITY=3.2100
    PRIMER_LEFT_4_PENALTY=0.721992
    PRIMER_RIGHT_4_PENALTY=0.050560
    PRIMER_LEFT_4_SEQUENCE=CAtttACACATTATTCCATTTCTtTTAtTC
    PRIMER_RIGHT_4_SEQUENCE=aTAGATTTTaAGAGTGGAaATAATTaTGAG
    PRIMER_LEFT_4=5294,30
    PRIMER_RIGHT_4=6180,30
    PRIMER_LEFT_4_TM=55.722
    PRIMER_RIGHT_4_TM=55.051
    PRIMER_LEFT_4_GC_PERCENT=23.333
    PRIMER_RIGHT_4_GC_PERCENT=23.333
    PRIMER_LEFT_4_SELF_ANY_TH=0.00
    PRIMER_RIGHT_4_SELF_ANY_TH=0.00
    PRIMER_LEFT_4_SELF_END_TH=0.00
    PRIMER_RIGHT_4_SELF_END_TH=0.00
    PRIMER_LEFT_4_HAIRPIN_TH=0.00
    PRIMER_RIGHT_4_HAIRPIN_TH=0.00
    PRIMER_LEFT_4_END_STABILITY=1.7500
    PRIMER_RIGHT_4_END_STABILITY=2.9000
    PRIMER_LEFT_5_PENALTY=0.543198
    PRIMER_RIGHT_5_PENALTY=3.641468
    PRIMER_LEFT_5_SEQUENCE=ATATATTTCAtATTTTTCAAAtTAtCCACC
    PRIMER_RIGHT_5_SEQUENCE=GAaTGAAaTGaAGaTaaAGAAGAGaAAGG
    PRIMER_LEFT_5=5764,30
    PRIMER_RIGHT_5=6701,29
    PRIMER_LEFT_5_TM=54.457
    PRIMER_RIGHT_5_TM=57.641
    PRIMER_LEFT_5_GC_PERCENT=20.000
    PRIMER_RIGHT_5_GC_PERCENT=31.034
    PRIMER_LEFT_5_SELF_ANY_TH=0.00
    PRIMER_RIGHT_5_SELF_ANY_TH=0.00
    PRIMER_LEFT_5_SELF_END_TH=0.00
    PRIMER_RIGHT_5_SELF_END_TH=0.00
    PRIMER_LEFT_5_HAIRPIN_TH=0.00
    PRIMER_RIGHT_5_HAIRPIN_TH=0.00
    PRIMER_LEFT_5_END_STABILITY=4.6100
    PRIMER_RIGHT_5_END_STABILITY=3.1100
    PRIMER_LEFT_6_PENALTY=0.195406
    PRIMER_RIGHT_6_PENALTY=2.701604
    PRIMER_LEFT_6_SEQUENCE=TCATCCCAtAAtAATTACATTAAATTtCAT
    PRIMER_RIGHT_6_SEQUENCE=aAaAaTGTaTGTGTGGAaaATTTTAGTTT
    PRIMER_LEFT_6=6269,30
    PRIMER_RIGHT_6=7171,29
    PRIMER_LEFT_6_TM=55.195
    PRIMER_RIGHT_6_TM=56.702
    PRIMER_LEFT_6_GC_PERCENT=20.000
    PRIMER_RIGHT_6_GC_PERCENT=20.690
    PRIMER_LEFT_6_SELF_ANY_TH=0.00
    PRIMER_RIGHT_6_SELF_ANY_TH=2.75
    PRIMER_LEFT_6_SELF_END_TH=0.00
    PRIMER_RIGHT_6_SELF_END_TH=0.00
    PRIMER_LEFT_6_HAIRPIN_TH=31.63
    PRIMER_RIGHT_6_HAIRPIN_TH=0.00
    PRIMER_LEFT_6_END_STABILITY=2.5700
    PRIMER_RIGHT_6_END_STABILITY=2.6600
    PRIMER_LEFT_7_PENALTY=0.018465
    PRIMER_LEFT_7_SEQUENCE=CAtCTTAtACCAtATtTCAttttTCTTtTt
    PRIMER_LEFT_7=6755,30
    PRIMER_LEFT_7_TM=54.982
    PRIMER_LEFT_7_GC_PERCENT=20.000
    PRIMER_LEFT_7_SELF_ANY_TH=0.00
    PRIMER_LEFT_7_SELF_END_TH=0.00
    PRIMER_LEFT_7_HAIRPIN_TH=0.00
    PRIMER_LEFT_7_END_STABILITY=1.9400
    =

Looks like we have two primers we can use for sequencing. Neither
overlaps a potential `CpG`, and both are outside of the C1q gene
sequence.

Next, we need to test them against the entire bisulfite-converted genome
to assess specificity.

# 6 GENOME BISULFITE CONVERSION

``` bash
# Load bash variables into memory
source .bashvars

tr '[:lower:]' '[:upper:]' < "${data_dir}/${genome_fasta}" | sed '/^>/!s/G/T/g' > "${data_dir}/${bisulfite_genome_fasta}"

head "${data_dir}/${bisulfite_genome_fasta}"
```

    >NC_052307.1 SALVELINUS NAMAYCUSH ISOLATE SENECA CHROMOSOME 1, SANAMA_1.0, WHOLE GENOME SHOTGUN SEQUENCE
    TTTTTAATTAATTATTTTACAATTTTTTTAAAACATAATTACTTTTTATTTTTATTCCACAATATAATTTAAATTACAAA
    TCTTTACTCCTATTATATTCTATATATTACCACATTATTATTCATAAAACTTATCTTTCAAACAAAAATCAACTTTTAAT
    ATACATTTCACATCTCATATATATCAATTCCATACATTTTTTTTTCTACTTCCTCATTTTTTTTTTTCTACTTCCTCATT
    TTTTTTTTTCTACTTCCTCATTTTTTTTTTTCTACTTCCTCATTTTTTTTTTTCTACTTCCTCATTTTTTTTTTTCTACT
    TCCTCATTTTTTTTTTTCTACTTCCTCATTTTTTTTTTTCTACTTCCTCATTTTCTACTTCCTCATTTTTTTTTTTCTAC
    ATCCTCATTTTTTTTTTTCTACTTCTTTTTTTTCTACTTCCTCATTTTTTTTTTTCTACATCCTCATTTTTTTTTTTCTA
    CTTCCTCATTTTTTTTTTTCTACTTCATCATTTTTTTTTTTCTACTTCTTTTTTTTCTACATCCTCATTTTTTTTTTTCT
    ACTTCCTCATTTTTTTTTTTCTACTTCCTCATTTTTTTTTTTCTACTTCCTCATTTTTTTTTTTCTACTTCCTCATTTTT
    TTTTTTCTACTTCTTTTTTTTCTACTTCCTCATTTTTTTTTTTCTACTTCCTCATTTTTTTTTTTCTACTTCTTTTTTTT

# 7 SPLIT GENOME

For some reason, I feel like the EMBOSS PrimerSearch tool will *not*
work on mulit-FastA files (or, it makes the results more difficult to
decipher?), but I’m just going off of memory and repeating a process
I’ve done previously.

This will split the bisulfite-converted genome multi-FastA file into
separate FastA files.

## 7.1 Split mulit-FastA file in to individual FastA files with PyFaidx

``` bash
# Load bash variables into memory
source .bashvars

# Make directory if it doesn't exist
mkdir --parents ${genome_fasta_splits_dir}

cd ${genome_fasta_splits_dir}

# Count sequences in FastA
echo "-------------------------------------------------------------------"
echo "NUMBER OF SEQUENCES IN ORIGINAL FASTA"
grep -c "^>" ../${bisulfite_genome_fasta}
echo "-------------------------------------------------------------------"
echo ""
echo ""



# Split FastA
${pyfaidx} \
--split-files \
../${bisulfite_genome_fasta}

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

# 8 PRIMER SEARCH WITH [EMOBSS PRIMERSEARCH](https://emboss.sourceforge.net/apps/cvs/emboss/apps/primersearch.html)

This will run [EMOBSS
PrimerSearch](https://emboss.sourceforge.net/apps/cvs/emboss/apps/primersearch.html)
against the genome to assess primer specificity.

## 8.1 Create EMOBSS PrimerSearch Primers File

Create a tab-delimited file to use with EMBOSS PrimerSearch.

``` bash
# Load bash variables into memory
source .bashvars

# Get primer info from Primer3 default format output
seq_id=$(grep "SEQUENCE_ID=" "${output_top}/primer3-t-primers-default-format.txt" | sed 's/SEQUENCE_ID=//')

left_primer=$(grep "PRIMER_LEFT_0_SEQUENCE=" "${output_top}/primer3-t-primers-default-format.txt" | sed 's/PRIMER_LEFT_0_SEQUENCE=//' | tr '[:lower:]' '[:upper:]')

right_primer=$(grep "PRIMER_RIGHT_6_SEQUENCE=" "${output_top}/primer3-t-primers-default-format.txt" | sed 's/PRIMER_RIGHT_6_SEQUENCE=//' | tr '[:lower:]' '[:upper:]')

# Create EMBOSS primer file
printf "%s\t%s\t%s\t\n" "${seq_id}" "${left_primer}" "${right_primer}" | tee "${output_top}/emboss-primers.txt"
```

    LOC120027825    TTATATTTCTTAACTACTAACCCTTATACA  AAAAATGTATGTGTGGAAAATTTTAGTTT   

## 8.2 Run EMBOSS PrimerSearch

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

    real    5m28.385s
    user    4m32.584s
    sys 1m14.440s

## 8.3 Check primer matches

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
        
        TTATATTTCTTAACTACTAACCCTTATACA hits forward strand at 19244321 with 0 mismatches
        AAAAATGTATGTGTGGAAAATTTTAGTTT hits reverse strand at [18091692] with 0 mismatches
        Amplimer length: 3907 bp
    ----------------------------------

# 9 SUMMARY

EMBOSS PrimerSearch has identified only a single location in the genome
where the forward and reverse primers match, as well as in the expected
location on chromosome `NC_052339`.

# 10 CITATIONS

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
