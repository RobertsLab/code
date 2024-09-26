20240202-e5-coral-sRNAseq-merging-explorations
================
Sam White
2024-02-02

- <a href="#1-intro" id="toc-1-intro">1 INTRO</a>
- <a href="#2-setup" id="toc-2-setup">2 SETUP</a>
  - <a href="#21-create-a-bash-variables-file"
    id="toc-21-create-a-bash-variables-file">2.1 Create a Bash variables
    file</a>
  - <a href="#22-download-raw-srnaseq-reads"
    id="toc-22-download-raw-srnaseq-reads">2.2 Download raw sRNAseq
    reads</a>
    - <a href="#221-verify-raw-read-checksums"
      id="toc-221-verify-raw-read-checksums">2.2.1 Verify raw read
      checksums</a>
  - <a href="#23-create-adapters-fasta"
    id="toc-23-create-adapters-fasta">2.3 Create adapters FastA</a>
  - <a href="#24-fastp" id="toc-24-fastp">2.4 <code>fastp</code></a>
    - <a href="#241-paired-end-trimming"
      id="toc-241-paired-end-trimming">2.4.1 Paired-end Trimming</a>
    - <a href="#242-single-end-trimming---r1-reads-only"
      id="toc-242-single-end-trimming---r1-reads-only">2.4.2 Single-end
      Trimming - R1 Reads Only</a>
  - <a href="#25-bbmerge" id="toc-25-bbmerge">2.5 BBMerge</a>
    - <a href="#251-raw-reads" id="toc-251-raw-reads">2.5.1 Raw reads</a>
    - <a href="#252-trimmed-polyg-only" id="toc-252-trimmed-polyg-only">2.5.2
      Trimmed PolyG Only</a>
    - <a href="#253-trimmed-adapters--polyg"
      id="toc-253-trimmed-adapters--polyg">2.5.3 Trimmed Adapters &amp;
      PolyG</a>
    - <a href="#254-trimmed-adapters--polyg-max-length-31bp"
      id="toc-254-trimmed-adapters--polyg-max-length-31bp">2.5.4 Trimmed
      Adapters &amp; PolyG Max Length 31bp</a>
  - <a href="#26-fastqc" id="toc-26-fastqc">2.6 FastQC</a>
  - <a href="#27-multiqc" id="toc-27-multiqc">2.7 MultiQC</a>

# 1 INTRO

Our previous sRNA-seq data analysis for this project revealed that the
R2 reads didn’t seem to match with anything (specifically, no clade
assignments when using
[miRTrace](https://github.com/friedlanderlab/mirtrace)), which seemed
odd. Since the insert sizes should be small for these libraries and the
read lengths kind of long (150bp PE), we’d expect the resulting reads to
overlap. And, in turn, R2 reads should be reverse complements of R1
reads (and vice versa). After discussing things with Azenta (company
which performed library construction and sequencing), they indicated
that their in-house pipeline discards the R2 reads because they come
from the opposite strand and most miRNAs are found on the R1 strand.

It seemed like a bit of waste to discard all the data that exists in the
R2 reads, so I wanted to explore merging the reads to see if we could
preserve the R2 data. Additionally, I came across
[XICRA](https://github.com/HCGB-IGTP/XICRA) (GitHub repo)
\[@sanchezherrero2021\] which describes using paired-end sequence data
to identify miRNA isoforms (isomiRs), suggesting that we should try to
retain the R2 reads to improve our miRNA analyses.

I decided to explore merging/trimming using two different sets of
software:

- [BBMerge](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)
  \[@bushnell2017\]

- [`fastp`](https://github.com/OpenGene/fastp)\[@chen2023; @chen2018\]

[BBMerge](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)
seemed like a good tool, as it’s specifically designed to handle merging
of Illumina PE sequencing reads. Plus, it’s part of the [BBTols
suite](https://sourceforge.net/projects/bbmap/), which has a number of
other tools, and packages like this usually have nice
cross-compatibility between the different tools.

[`fastp`](https://github.com/OpenGene/fastp) has both trimming *and*
merging built into the tools itself, which I hadn’t previously realized
(but, had never explored, since this isn’t how we usually work with PE
sequencing data). I was simply going to use
[`fastp`](https://github.com/OpenGene/fastp) for trimming, as it’s the
trimmer thta I’m most familiar with. I honestly thought I’d trim with
[`fastp`](https://github.com/OpenGene/fastp), and then merge with
[BBMerge](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/),
but the fact that [`fastp`](https://github.com/OpenGene/fastp) had it
built-in could make my life easier.

So, with all of that out of the way, let’s start the analysis. This will
only use a single set of PE sequencing reads, as an initial exploration.
We’ll decide on the best approach and then use that approach to run all
the sequencing data through the chosen “pipeline.”

# 2 SETUP

## 2.1 Create a Bash variables file

This allows usage of Bash variables across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "## Data URL"
echo 'export raw_reads_url="https://owl.fish.washington.edu/nightingales/P_evermanni/30-852430235/"'

echo "# Data directories"
echo 'export fastqc_out_dir="./fastqc"'
echo 'export raw_reads_dir="./raw"'
echo 'export trimmed_fastqs_dir="./trimmed-reads"'
echo ""

echo "## NEB nebnext-small-rna-library-prep-set-for-illumina adapters"
echo 'export first_adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"'
echo 'export second_adapter="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"'
echo ""

echo "# Input/output files"
echo 'export NEB_adapters_fasta=NEB-adapters.fasta'
echo ""

echo "# Paths to programs"
echo 'export programs_dir="/home/shared"'
echo 'export bbmerge="${programs_dir}/bbmap-39.06/bbmerge.sh"'
echo 'export fastp="${programs_dir}/fastp"'
echo 'export fastqc="${programs_dir}/FastQC-0.12.1/fastqc"'
echo 'export multiqc=/home/sam/programs/mambaforge/bin/multiqc'
echo ""

echo "# Set number of CPUs to use"
echo 'export threads=46'
echo ""

echo "# Initialize arrays"
echo 'export trimmed_fastqs_array=()'


} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    ## Data URL
    export raw_reads_url="https://owl.fish.washington.edu/nightingales/P_evermanni/30-852430235/"
    # Data directories
    export fastqc_out_dir="./fastqc"
    export raw_reads_dir="./raw"
    export trimmed_fastqs_dir="./trimmed-reads"

    ## NEB nebnext-small-rna-library-prep-set-for-illumina adapters
    export first_adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    export second_adapter="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"

    # Input/output files
    export NEB_adapters_fasta=NEB-adapters.fasta

    # Paths to programs
    export programs_dir="/home/shared"
    export bbmerge="${programs_dir}/bbmap-39.06/bbmerge.sh"
    export fastp="${programs_dir}/fastp"
    export fastqc="${programs_dir}/FastQC-0.12.1/fastqc"
    export multiqc=/home/sam/programs/mambaforge/bin/multiqc

    # Set number of CPUs to use
    export threads=46

    # Initialize arrays
    export trimmed_fastqs_array=()

## 2.2 Download raw sRNAseq reads

Reads are downloaded from
<https://owl.fish.washington.edu/nightingales/P_evermanni/30-852430235/>

The `--cut-dirs 3` command cuts the preceding directory structure
(i.e. `nightingales/P_evermanni/30-852430235/`) so that we just end up
with the reads.

``` bash
# Load bash variables into memory
source .bashvars

mkdir --parents "${raw_reads_dir}" "${trimmed_fastqs_dir}"

wget \
--directory-prefix ${raw_reads_dir} \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 3 \
--no-host-directories \
--no-parent \
--quiet \
--accept "sRNA-POR-79*,checksums.md5" ${raw_reads_url}

ls -lh "${raw_reads_dir}"
```

    total 1.8G
    -rw-r--r-- 1 sam sam  798 May 17  2023 checksums.md5
    -rw-r--r-- 1 sam sam 899M May 17  2023 sRNA-POR-79-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 sam sam 916M May 17  2023 sRNA-POR-79-S1-TP2_R2_001.fastq.gz

### 2.2.1 Verify raw read checksums

``` bash
# Load bash variables into memory
source .bashvars

cd "${raw_reads_dir}"

# Checksums file contains other files, so this just looks for the sRNAseq files.
grep "sRNA-POR-79" checksums.md5 | md5sum --check
```

    sRNA-POR-79-S1-TP2_R1_001.fastq.gz: OK
    sRNA-POR-79-S1-TP2_R2_001.fastq.gz: OK

## 2.3 Create adapters FastA

``` bash
# Load bash variables into memory
source .bashvars

echo "Creating adapters FastA."
echo ""
adapter_count=0

# Check for adapters file first
# Then create adapters file if doesn't exist
if [ -f "${NEB_adapters_fasta}" ]; then
  echo "${NEB_adapters_fasta} already exists. Nothing to do."
else
  for adapter in "${first_adapter}" "${second_adapter}"
  do
    adapter_count=$((adapter_count + 1))
    printf ">%s\n%s\n" "adapter_${adapter_count}" "${adapter}"
  done >> "${NEB_adapters_fasta}"
fi

echo ""
echo "Adapters FastA:"
echo ""
cat "${NEB_adapters_fasta}"
echo ""
```

    Creating adapters FastA.


    Adapters FastA:

    >adapter_1
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    >adapter_2
    GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT

## 2.4 [`fastp`](https://github.com/OpenGene/fastp)

[`fastp`](https://github.com/OpenGene/fastp) has many built-in
functions, including before/after trimming plots, and read 1 and read 2
merging.

<div class="callout-note">

Setting `--overlap_len_require 17` forces
[`fastp`](https://github.com/OpenGene/fastp) to look at reads starting
at base 17, instead of the default 30 and results in more accurate
analyses, since these are sRNA-seq libraries where the expected/desired
read lengths will be 17 - 30bp in length.

</div>

### 2.4.1 Paired-end Trimming

The paired end trimming will examine how some of the trimming options
impact sequencing read lengths, quality, and subsequent merging.

#### 2.4.1.1 PolyG Trimming Only

``` bash
# Load bash variables into memory
source .bashvars

"${fastp}" \
--in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz \
--in2 ./raw/sRNA-POR-79-S1-TP2_R2_001.fastq.gz \
--out1 ./trimmed-reads/fastp-polyG_only.R1.fq.gz \
--out2 ./trimmed-reads/fastp-polyG_only.R2.fq.gz \
--disable_adapter_trimming \
--trim_poly_g \
--overlap_len_require 17 \
--thread ${threads} \
--html "fastp-polyG-only.html" \
--json "fastp-polyG-only.json" \
--report_title "fastp-polyG-only"
```

    WARNING: fastp uses up to 16 threads although you specified 46
    Read1 before filtering:
    total reads: 19762559
    total bases: 2964383850
    Q20 bases: 2604771220(87.8689%)
    Q30 bases: 2427856243(81.9009%)

    Read2 before filtering:
    total reads: 19762559
    total bases: 2964383850
    Q20 bases: 2619314070(88.3595%)
    Q30 bases: 2408732242(81.2557%)

    Read1 after filtering:
    total reads: 19450041
    total bases: 1698288792
    Q20 bases: 1478970498(87.0859%)
    Q30 bases: 1370182460(80.6802%)

    Read2 after filtering:
    total reads: 19450041
    total bases: 1849068926
    Q20 bases: 1690862447(91.444%)
    Q30 bases: 1596965342(86.3659%)

    Filtering result:
    reads passed filter: 38900082
    reads failed due to low quality: 595592
    reads failed due to too many N: 40
    reads failed due to too short: 29404

    Duplication rate: 10.4184%

    Insert size peak (evaluated by paired-end reads): 29

    JSON report: fastp-polyG-only.json
    HTML report: fastp-polyG-only.html

    /home/shared/fastp --in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz --in2 ./raw/sRNA-POR-79-S1-TP2_R2_001.fastq.gz --out1 ./trimmed-reads/fastp-polyG_only.R1.fq.gz --out2 ./trimmed-reads/fastp-polyG_only.R2.fq.gz --disable_adapter_trimming --trim_poly_g --overlap_len_require 17 --thread 46 --html fastp-polyG-only.html --json fastp-polyG-only.json --report_title fastp-polyG-only 
    fastp v0.23.2, time used: 43 seconds

#### 2.4.1.2 Adapters & PolyG Trimming

``` bash
# Load bash variables into memory
source .bashvars

"${fastp}" \
--in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz \
--in2 ./raw/sRNA-POR-79-S1-TP2_R2_001.fastq.gz \
--out1 ./trimmed-reads/fastp-adapters-polyG.R1.fq.gz \
--out2 ./trimmed-reads/fastp-adapters-polyG.R2.fq.gz \
--adapter_fasta NEB-adapters.fasta \
--trim_poly_g \
--overlap_len_require 17 \
--thread ${threads} \
--html "fastp-adapters-polyG.html" \
--json "fastp-adapters-polyG.json" \
--report_title "fastp-adapters-polyG-only"
```

    WARNING: fastp uses up to 16 threads although you specified 46
    Read1 before filtering:
    total reads: 19762559
    total bases: 2964383850
    Q20 bases: 2604771220(87.8689%)
    Q30 bases: 2427856243(81.9009%)

    Read2 before filtering:
    total reads: 19762559
    total bases: 2964383850
    Q20 bases: 2619314070(88.3595%)
    Q30 bases: 2408732242(81.2557%)

    Read1 after filtering:
    total reads: 16734573
    total bases: 489646406
    Q20 bases: 480139354(98.0584%)
    Q30 bases: 466862298(95.3468%)

    Read2 after filtering:
    total reads: 16734573
    total bases: 492483642
    Q20 bases: 484765577(98.4328%)
    Q30 bases: 469308001(95.2941%)

    Filtering result:
    reads passed filter: 33469146
    reads failed due to low quality: 232788
    reads failed due to too many N: 0
    reads failed due to too short: 5823184
    reads with adapter trimmed: 30461598
    bases trimmed due to adapters: 2013836548

    Duplication rate: 10.4184%

    Insert size peak (evaluated by paired-end reads): 29

    JSON report: fastp-adapters-polyG.json
    HTML report: fastp-adapters-polyG.html

    /home/shared/fastp --in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz --in2 ./raw/sRNA-POR-79-S1-TP2_R2_001.fastq.gz --out1 ./trimmed-reads/fastp-adapters-polyG.R1.fq.gz --out2 ./trimmed-reads/fastp-adapters-polyG.R2.fq.gz --adapter_fasta NEB-adapters.fasta --trim_poly_g --overlap_len_require 17 --thread 46 --html fastp-adapters-polyG.html --json fastp-adapters-polyG.json --report_title fastp-adapters-polyG-only 
    fastp v0.23.2, time used: 32 seconds

#### 2.4.1.3 PolyG Trimming Max Length 31

This can’t be run because we know from the polyG only trimming, that
most read lengths are \>89bp. As such, it will end up discarding all the
reads…

#### 2.4.1.4 Adapters & PolyG Trimming Max Length 31

<div class="callout-note">

The max length is based on the `fastp` insert peak size from the adapter
and polyG trimming results, and previous evaluation of mean read lengths
via
[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
and [`MultiQC`](https://multiqc.info/).

</div>

``` bash
# Load bash variables into memory
source .bashvars

"${fastp}" \
--in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz \
--in2 ./raw/sRNA-POR-79-S1-TP2_R2_001.fastq.gz \
--out1 ./trimmed-reads/fastp-adapters-polyG-31bp.R1.fq.gz \
--out2 ./trimmed-reads/fastp-adapters-polyG-31bp.R2.fq.gz \
--adapter_fasta NEB-adapters.fasta \
--trim_poly_g \
--overlap_len_require 17 \
--length_limit 31 \
--thread ${threads} \
--html "fastp-adapters-polyG-31bp.html" \
--json "fastp-adapters-polyG-31bp.json" \
--report_title "fastp-adapters-polyG-31bp"
```

    WARNING: fastp uses up to 16 threads although you specified 46
    Read1 before filtering:
    total reads: 19762559
    total bases: 2964383850
    Q20 bases: 2604771220(87.8689%)
    Q30 bases: 2427856243(81.9009%)

    Read2 before filtering:
    total reads: 19762559
    total bases: 2964383850
    Q20 bases: 2619314070(88.3595%)
    Q30 bases: 2408732242(81.2557%)

    Read1 after filtering:
    total reads: 13190855
    total bases: 339210880
    Q20 bases: 335529961(98.9149%)
    Q30 bases: 327446980(96.532%)

    Read2 after filtering:
    total reads: 13190855
    total bases: 338558013
    Q20 bases: 335504115(99.098%)
    Q30 bases: 325352134(96.0994%)

    Filtering result:
    reads passed filter: 26381710
    reads failed due to low quality: 232788
    reads failed due to too many N: 0
    reads failed due to too short: 5354066
    reads failed due to too long: 7556554
    reads with adapter trimmed: 30461598
    bases trimmed due to adapters: 2013836548

    Duplication rate: 10.4184%

    Insert size peak (evaluated by paired-end reads): 29

    JSON report: fastp-adapters-polyG-31bp.json
    HTML report: fastp-adapters-polyG-31bp.html

    /home/shared/fastp --in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz --in2 ./raw/sRNA-POR-79-S1-TP2_R2_001.fastq.gz --out1 ./trimmed-reads/fastp-adapters-polyG-31bp.R1.fq.gz --out2 ./trimmed-reads/fastp-adapters-polyG-31bp.R2.fq.gz --adapter_fasta NEB-adapters.fasta --trim_poly_g --overlap_len_require 17 --length_limit 31 --thread 46 --html fastp-adapters-polyG-31bp.html --json fastp-adapters-polyG-31bp.json --report_title fastp-adapters-polyG-31bp 
    fastp v0.23.2, time used: 31 seconds

#### 2.4.1.5 Adapters & PolyG Trimming Max Length 31 With Merge

This will attempt to merge the R1 and R2 reads.

``` bash
# Load bash variables into memory
source .bashvars

"${fastp}" \
--in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz \
--in2 ./raw/sRNA-POR-79-S1-TP2_R2_001.fastq.gz \
--adapter_fasta NEB-adapters.fasta \
--trim_poly_g \
--overlap_len_require 17 \
--length_limit 31 \
--merge \
--merged_out ./trimmed-reads/fastp-adapters-polyG-31bp-merged.fq.gz \
--thread ${threads} \
--html "fastp-adapters-polyG-31bp-merged.html" \
--json "fastp-adapters-polyG-31bp-merged.json" \
--report_title "fastp-adapters-polyG-31bp-merged"
```

    WARNING: fastp uses up to 16 threads although you specified 46
    Read1 before filtering:
    total reads: 19762559
    total bases: 2964383850
    Q20 bases: 2604771220(87.8689%)
    Q30 bases: 2427856243(81.9009%)

    Read2 before filtering:
    total reads: 19762559
    total bases: 2964383850
    Q20 bases: 2619314070(88.3595%)
    Q30 bases: 2408732242(81.2557%)

    Merged and filtered:
    total reads: 12259576
    total bases: 324356716
    Q20 bases: 320917135(98.9396%)
    Q30 bases: 313148708(96.5445%)

    Filtering result:
    reads passed filter: 26506174
    reads failed due to low quality: 179320
    reads failed due to too many N: 0
    reads failed due to too short: 5774946
    reads failed due to too long: 7064678
    reads with adapter trimmed: 30461598
    bases trimmed due to adapters: 2013836548
    reads corrected by overlap analysis: 484376
    bases corrected by overlap analysis: 766038

    Duplication rate: 10.4184%

    Insert size peak (evaluated by paired-end reads): 29

    Read pairs merged: 12259576
    % of original read pairs: 62.0344%
    % in reads after filtering: 100%


    JSON report: fastp-adapters-polyG-31bp-merged.json
    HTML report: fastp-adapters-polyG-31bp-merged.html

    /home/shared/fastp --in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz --in2 ./raw/sRNA-POR-79-S1-TP2_R2_001.fastq.gz --adapter_fasta NEB-adapters.fasta --trim_poly_g --overlap_len_require 17 --length_limit 31 --merge --merged_out ./trimmed-reads/fastp-adapters-polyG-31bp-merged.fq.gz --thread 46 --html fastp-adapters-polyG-31bp-merged.html --json fastp-adapters-polyG-31bp-merged.json --report_title fastp-adapters-polyG-31bp-merged 
    fastp v0.23.2, time used: 43 seconds

### 2.4.2 Single-end Trimming - R1 Reads Only

Using just the R1 reads, per Azenta’s recommendation.

#### 2.4.2.1 Auto-adapters & PolyG Trimming

``` bash
# Load bash variables into memory
source .bashvars

"${fastp}" \
--in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz \
--out1 ./trimmed-reads/fastp-R1-auto_adapters-polyG.fq.gz \
--trim_poly_g \
--thread ${threads} \
--html "fastp-R1-auto_adapters-polyG.html" \
--json "fastp-R1-auto_adapters-polyG.json" \
--report_title "fastp-R1-auto_adapters-polyG"
```

    Detecting adapter sequence for read1...
    >Illumina TruSeq Adapter Read 1
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

    WARNING: fastp uses up to 16 threads although you specified 46
    Read1 before filtering:
    total reads: 19762559
    total bases: 2964383850
    Q20 bases: 2604771220(87.8689%)
    Q30 bases: 2427856243(81.9009%)

    Read1 after filtering:
    total reads: 16984592
    total bases: 540469705
    Q20 bases: 522006375(96.5838%)
    Q30 bases: 504657126(93.3738%)

    Filtering result:
    reads passed filter: 16984592
    reads failed due to low quality: 30923
    reads failed due to too many N: 1
    reads failed due to too short: 2747043
    reads with adapter trimmed: 18545632
    bases trimmed due to adapters: 1166285709

    Duplication rate (may be overestimated since this is SE data): 29.621%

    JSON report: fastp-R1-auto_adapters-polyG.json
    HTML report: fastp-R1-auto_adapters-polyG.html

    /home/shared/fastp --in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz --out1 ./trimmed-reads/fastp-R1-auto_adapters-polyG.fq.gz --trim_poly_g --thread 46 --html fastp-R1-auto_adapters-polyG.html --json fastp-R1-auto_adapters-polyG.json --report_title fastp-R1-auto_adapters-polyG 
    fastp v0.23.2, time used: 31 seconds

#### 2.4.2.2 Auto-adapters & PolyG Trimming Max Length 31

<div class="callout-note">

The max length is based on the `fastp` insert peak size from the adapter
and polyG trimming results, and previous evaluation of mean read lengths
via
[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
and [`MultiQC`](https://multiqc.info/).

</div>

``` bash
# Load bash variables into memory
source .bashvars

"${fastp}" \
--in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz \
--out1 ./trimmed-reads/fastp-R1-31bp-auto_adapters-polyG.fq.gz \
--trim_poly_g \
--length_limit 31 \
--thread ${threads} \
--html "fastp-R1-auto_adapters-polyG-31bp.html" \
--json "fastp-R1-auto_adapters-polyG-31bp.json" \
--report_title "fastp-R1-auto_adapters-polyG-31bp"
```

    Detecting adapter sequence for read1...
    >Illumina TruSeq Adapter Read 1
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

    WARNING: fastp uses up to 16 threads although you specified 46
    Read1 before filtering:
    total reads: 19762559
    total bases: 2964383850
    Q20 bases: 2604771220(87.8689%)
    Q30 bases: 2427856243(81.9009%)

    Read1 after filtering:
    total reads: 12801426
    total bases: 328647991
    Q20 bases: 325069951(98.9113%)
    Q30 bases: 317201339(96.517%)

    Filtering result:
    reads passed filter: 12801426
    reads failed due to low quality: 30923
    reads failed due to too many N: 1
    reads failed due to too short: 2747043
    reads failed due to too long: 4183166
    reads with adapter trimmed: 18545632
    bases trimmed due to adapters: 1166285709

    Duplication rate (may be overestimated since this is SE data): 29.621%

    JSON report: fastp-R1-auto_adapters-polyG-31bp.json
    HTML report: fastp-R1-auto_adapters-polyG-31bp.html

    /home/shared/fastp --in1 ./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz --out1 ./trimmed-reads/fastp-R1-31bp-auto_adapters-polyG.fq.gz --trim_poly_g --length_limit 31 --thread 46 --html fastp-R1-auto_adapters-polyG-31bp.html --json fastp-R1-auto_adapters-polyG-31bp.json --report_title fastp-R1-auto_adapters-polyG-31bp 
    fastp v0.23.2, time used: 31 seconds

## 2.5 BBMerge

### 2.5.1 Raw reads

[BBMerge](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)
documentations actually indicates that merging performs best using
*raw*, untrimmed reads, so we’ll try that out first…

``` bash
# Load bash variables into memory
source .bashvars
${bbmerge} \
in1=./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz \
in2=./raw/sRNA-POR-79-S1-TP2_R2_001.fastq.gz \
adapters=./NEB-adapters.fasta \
mininsert=17 \
out=./trimmed-reads/sRNA-POR-79-S1-TP2-raw-bbmerge.fq.gz
```

    java -ea -Xmx1000m -Xms1000m -Djava.library.path=/home/shared/bbmap-39.06/jni/ -cp /home/shared/bbmap-39.06/current/ jgi.BBMerge in1=./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz in2=./raw/sRNA-POR-79-S1-TP2_R2_001.fastq.gz adapters=./NEB-adapters.fasta mininsert=17 out=./trimmed-reads/sRNA-POR-79-S1-TP2-raw-bbmerge.fq.gz
    Executing jgi.BBMerge [in1=./raw/sRNA-POR-79-S1-TP2_R1_001.fastq.gz, in2=./raw/sRNA-POR-79-S1-TP2_R2_001.fastq.gz, adapters=./NEB-adapters.fasta, mininsert=17, out=./trimmed-reads/sRNA-POR-79-S1-TP2-raw-bbmerge.fq.gz]
    Version 39.06

    Writing mergable reads merged.
    Started output threads.
    Total time: 56.412 seconds.

    Pairs:                  19762559
    Joined:                 151908      0.769%
    Ambiguous:              0           0.000%
    No Solution:            19610645    99.231%
    Too Short:              6           0.000%
    Adapters Expected:      60129       0.152%
    Adapters Found:         58428       0.148%

    Avg Insert:             131.5
    Standard Deviation:     62.2
    Mode:                   156

    Insert range:           17 - 289
    90th percentile:        203
    75th percentile:        156
    50th percentile:        156
    25th percentile:        72
    10th percentile:        32

### 2.5.2 Trimmed PolyG Only

Running
[BBMerge](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)
on just raw reads yielded poor results (\>99% of reads could *not* be
merged). Additionally, the average insert size was determined to be
131.5bp, which is not going to be useful for sRNA-seq.

I know that there’s a very large stretch of polyG sequence in the reads
(likely due to the small insert size and “long” read length - Gs get
added as a “placeholder” base each cycle that exceeds the insert size),
so let’s try merging reads *after* polyG sequence has been trimmed.

``` bash
# Load bash variables into memory
source .bashvars

${bbmerge} \
in1=./trimmed-reads/fastp-polyG_only.R1.fq.gz \
in2=./trimmed-reads/fastp-polyG_only.R2.fq.gz \
adapters=./NEB-adapters.fasta \
mininsert=17 \
out=./trimmed-reads/sRNA-POR-79-S1-TP2-fastp-polyG_only-bbmerge.fq.gz
```

    java -ea -Xmx1000m -Xms1000m -Djava.library.path=/home/shared/bbmap-39.06/jni/ -cp /home/shared/bbmap-39.06/current/ jgi.BBMerge in1=./trimmed-reads/fastp-polyG_only.R1.fq.gz in2=./trimmed-reads/fastp-polyG_only.R2.fq.gz adapters=./NEB-adapters.fasta mininsert=17 out=./trimmed-reads/sRNA-POR-79-S1-TP2-fastp-polyG_only-bbmerge.fq.gz
    Executing jgi.BBMerge [in1=./trimmed-reads/fastp-polyG_only.R1.fq.gz, in2=./trimmed-reads/fastp-polyG_only.R2.fq.gz, adapters=./NEB-adapters.fasta, mininsert=17, out=./trimmed-reads/sRNA-POR-79-S1-TP2-fastp-polyG_only-bbmerge.fq.gz]
    Version 39.06

    Writing mergable reads merged.
    Started output threads.
    Total time: 52.363 seconds.

    Pairs:                  19450041
    Joined:                 12649853    65.038%
    Ambiguous:              0           0.000%
    No Solution:            5981377     30.753%
    Too Short:              818811      4.210%
    Adapters Expected:      12781441    32.857%
    Adapters Found:         12553911    32.272%

    Avg Insert:             30.0
    Standard Deviation:     14.6
    Mode:                   29

    Insert range:           17 - 289
    90th percentile:        37
    75th percentile:        31
    50th percentile:        29
    25th percentile:        25
    10th percentile:        21

The results from using polyG-trimmed reads are certainly much better
than just merging raw reads (12,649,853 reads could be joined this
time). Mean insert size is now determined to be 30bp, which is what we
might expect. However, the insert range (17 - 289) is still concerning,
as reads \>30bp likely aren’t usable. Adapter sequences are still
present, which are certainly contributing to some of the inability to
merge, as well as the large read lengths. So, let’s try merging reads
that have been trimmed of both polyG and adapter sequences.

### 2.5.3 Trimmed Adapters & PolyG

``` bash
# Load bash variables into memory
source .bashvars

${bbmerge} \
in1=./trimmed-reads/fastp-adapters-polyG.R1.fq.gz \
in2=./trimmed-reads/fastp-adapters-polyG.R2.fq.gz \
adapters=./NEB-adapters.fasta \
mininsert=17 \
out=./trimmed-reads/sRNA-POR-79-S1-TP2-fastp-adapters-polyG-bbmerge.fq.gz
```

    java -ea -Xmx1000m -Xms1000m -Djava.library.path=/home/shared/bbmap-39.06/jni/ -cp /home/shared/bbmap-39.06/current/ jgi.BBMerge in1=./trimmed-reads/fastp-adapters-polyG.R1.fq.gz in2=./trimmed-reads/fastp-adapters-polyG.R2.fq.gz adapters=./NEB-adapters.fasta mininsert=17 out=./trimmed-reads/sRNA-POR-79-S1-TP2-fastp-adapters-polyG-bbmerge.fq.gz
    Executing jgi.BBMerge [in1=./trimmed-reads/fastp-adapters-polyG.R1.fq.gz, in2=./trimmed-reads/fastp-adapters-polyG.R2.fq.gz, adapters=./NEB-adapters.fasta, mininsert=17, out=./trimmed-reads/sRNA-POR-79-S1-TP2-fastp-adapters-polyG-bbmerge.fq.gz]
    Version 39.06

    Writing mergable reads merged.
    Started output threads.
    Total time: 27.985 seconds.

    Pairs:                  16734573
    Joined:                 15533284    92.822%
    Ambiguous:              0           0.000%
    No Solution:            603776      3.608%
    Too Short:              597513      3.571%
    Adapters Expected:      28222       0.084%
    Adapters Found:         18307       0.055%

    Avg Insert:             29.4
    Standard Deviation:     13.2
    Mode:                   29

    Insert range:           17 - 289
    90th percentile:        36
    75th percentile:        31
    50th percentile:        29
    25th percentile:        25
    10th percentile:        21

Well, trimming adapters certainly made a difference, as well! Now,
15,533,284 reads get joined (as opposed to 12,649,853 when using just
polyG-trimmed reads). The average insert size changes very little, but,
oddly, the insert range remains the same (17 - 289).

For kicks, let’s throw the reads that have been trimmed to 31bp into
[BBMerge](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)
and see what happens.

### 2.5.4 Trimmed Adapters & PolyG Max Length 31bp

``` bash
# Load bash variables into memory
source .bashvars

${bbmerge} \
in1=./trimmed-reads/fastp-adapters-polyG-31bp.R1.fq.gz \
in2=./trimmed-reads/fastp-adapters-polyG-31bp.R2.fq.gz \
adapters=./NEB-adapters.fasta \
mininsert=17 \
out=./trimmed-reads/sRNA-POR-79-S1-TP2-fastp-adapters-polyG-31bp-bbmerge.fq.gz
```

    java -ea -Xmx1000m -Xms1000m -Djava.library.path=/home/shared/bbmap-39.06/jni/ -cp /home/shared/bbmap-39.06/current/ jgi.BBMerge in1=./trimmed-reads/fastp-adapters-polyG-31bp.R1.fq.gz in2=./trimmed-reads/fastp-adapters-polyG-31bp.R2.fq.gz adapters=./NEB-adapters.fasta mininsert=17 out=./trimmed-reads/sRNA-POR-79-S1-TP2-fastp-adapters-polyG-31bp-bbmerge.fq.gz
    Executing jgi.BBMerge [in1=./trimmed-reads/fastp-adapters-polyG-31bp.R1.fq.gz, in2=./trimmed-reads/fastp-adapters-polyG-31bp.R2.fq.gz, adapters=./NEB-adapters.fasta, mininsert=17, out=./trimmed-reads/sRNA-POR-79-S1-TP2-fastp-adapters-polyG-31bp-bbmerge.fq.gz]
    Version 39.06

    Writing mergable reads merged.
    Started output threads.
    Total time: 21.354 seconds.

    Pairs:                  13190855
    Joined:                 12346728    93.601%
    Ambiguous:              0           0.000%
    No Solution:            281891      2.137%
    Too Short:              562236      4.262%
    Adapters Expected:      4           0.000%
    Adapters Found:         4           0.000%

    Avg Insert:             26.3
    Standard Deviation:     4.0
    Mode:                   29

    Insert range:           17 - 35
    90th percentile:        31
    75th percentile:        29
    50th percentile:        28
    25th percentile:        24
    10th percentile:        20

This resulted in 12,346,728 reads being joined. This is on par with the
12,649,853 reads joined when using just polyG-trimmed reads.
Additionally, we see the average insert size is 26.3bp, which is even
more inline with what we’d expect for miRNAs (22-24bp), and the insert
range is also inline with what we’d expect for these libraries (17 -
35bp). However, this is about 20% fewer reads than the 15,533,284 reads
when using polyG- and adapter-trimmed reads.

Maybe I’ll test out the subsequent downstream analyses using these
different merging results…

## 2.6 FastQC

``` bash
# Load bash variables into memory
source .bashvars

# Populate array with FastQ files
fastq_array=(./trimmed-reads/*.fq.gz)

# Pass array contents to new variable
fastqc_list=$(echo "${fastq_array[*]}")

# Run FastQC
# NOTE: Do NOT quote ${fastqc_list}
${fastqc} \
--threads ${threads} \
--quiet \
--outdir ./ \
${fastqc_list}
```

    application/gzip
    application/gzip
    application/gzip
    application/gzip
    application/gzip
    application/gzip
    application/gzip
    application/gzip
    application/gzip
    application/gzip
    application/gzip
    application/gzip
    application/gzip

## 2.7 MultiQC

``` bash
# Load bash variables into memory
source .bashvars

${multiqc} .
```

      /// MultiQC 🔍 | v1.14

    |           multiqc | MultiQC Version v1.25 now available!
    |           multiqc | Search path : /home/shared/8TB_HDD_01/sam/gitrepos/RobertsLab/code/r_projects/sam/20240202-e5-coral-sRNAseq-merging-explorations/code
    |         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 57/57  
    |            fastqc | Found 13 reports
    |           multiqc | Compressing plot data
    |           multiqc | Report      : multiqc_report.html
    |           multiqc | Data        : multiqc_data
    |           multiqc | MultiQC complete

Check out the [`MultiQC`](https://multiqc.info/) results here (HTML):

- [MultiQC Report](./multiqc_report.html)
