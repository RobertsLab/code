#!/bin/bash
## Job Name - can be changed
#SBATCH --job-name=BISMARK
## Allocation Definition - confirm correctness
#SBATCH --account=coenv
#SBATCH --partition=coenv
## Resources
## Nodes (often you will only use 1)
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=30-00:00:00
## Memory per node
#SBATCH --mem=120G
## email notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER@uw.edu
## Specify the working directory for this job
#SBATCH --workdir=
# Exit script if a command fails
set -e

### USER NEEDS TO SET VARIABLES FOR THE FOLLOWING:
# Full path to directory with sequencing reads
reads_dir=""

# Full path to bisulftie-converted genome directory
genome_dir=""

# Enter y (for yes) or n (for no) between the quotes.
# Yes - Whole genome bisulfite sequencing
# No - Reduced genome bisulfite sequencing (e.g. RRBS, MBD)
deduplicated=""
####################################################
# DO NOT EDIT BELOW THIS LINE

# Evaluate user-edited variables to make sure they have been filled
[ ! -z ${deduplicated} ] \
&& echo "The deduplicated variable is not defined. Please edit the SBATCH script and add y or n to deduplicated variable."

[ ! -z ${genome_dir} ] \
&& echo "The bisulfite genome directory path has not been set. Please edit the SBATCH script."

[ ! -z ${reads_dir} ] \
&& echo "The reads directory paht has not been set. Please edit the SBATCH script."

# Directories and programs
wd=$(pwd)
bismark_dir="/gscratch/srlab/programs/Bismark-0.21.0"
bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/"
samtools="/gscratch/srlab/programs/samtools-1.9/samtools"
threads="28"
reads_list="input_fastqs.txt"

# Create list of input FastQ files for easier confirmation.
for fastq in ${reads_dir}/*.fq*
do
  echo ${fastq} >> ${reads_list}
done

# Initialize arrays
R1_array=()
R2_array=()
bam_array=()
dedup_bam_array=()

# Check for paired-end
# Capture grep output
# >0 means single-end reads
grep "_R2_" ${reads_list}
paired=$?

# Concatenate R1 reads and generate lists of FastQs
for fastq in ${reads_dir}/*R1*.gz
do
  cat ${fastq} >> ${R1}
done

## Save FastQ files to arrays
R1_array=(${reads_dir}/*_R1_*.fq*)


if [ ${paired} -eq 0 ]; then
  # Concatenate R2 reads
  R2_array=(${reads_dir}/*_R2_*.fq*)
  # Concatenate R2 reads and generate lists of FastQs
  for fastq in ${reads_dir}/*R2*.gz
    do
      cat ${fastq} >> ${R2}
  done
  # Run bismark using bisulftie-converted genome
  # Generates a set of BAM files as outputs
  # Records stderr to a file for easy viewing of Bismark summary info
  ${bismark_dir}/bismark \
  --path_to_bowtie ${bowtie2_dir} \
  --genome ${genome} \
  --non_directional \
  -p 28 \
  -1 ${R1} \
  -2 ${R2} \
  2> bismark_summary.txt
else
  ${bismark_dir}/bismark \
  --path_to_bowtie ${bowtie2_dir} \
  --genome ${genome} \
  --non_directional \
  -p 28 \
  ${R1} \
  2> bismark_summary.txt
fi

# Sort Bismark BAM files (failsafe for deduplication)
find *.bam \
| xargs basename -s .bam \
| xargs -I{} \
${samtools} sort \
--threads ${threads} \
-n {}.bam \
-o {}.sorted.bam

# Deduplication
find *sorted.bam \
| xargs basename -s .bam \
| xargs -I{} \
${bismark_dir}/deduplicate_bismark \
--paired \
--samtools_path=${samtools} \
{}.bam

# Methylation extraction
# Extracts methylation info from BAM files produced by Bismark
# Options to created a bedgraph file, counts, remove spaces from names
# and to use the "scaffolds" setting.
${bismark_dir}/bismark_methylation_extractor \
--bedGraph \
--cytosine_report \
--counts \
--scaffolds \
--remove_spaces \
--multicore ${threads} \
--buffer_size 75% \
*.bam


# Methylation extraction
# Extracts methylation info from deduplicated BAM files produced by Bismark
# Options to created a bedgraph file, counts, remove spaces from names
# and to use the "scaffolds" setting.
${bismark_dir}/bismark_methylation_extractor \
--bedGraph \
--cytosine_report \
--counts \
--scaffolds \
--remove_spaces \
--multicore ${threads} \
--buffer_size 75% \
*deduplicated.bam


# Bismark processing report
# Generates HTML reports from previously created files
${bismark_dir}/bismark2report

#Bismark summary report
# Generates HTML summary reports from previously created files
${bismark_dir}/bismark2summary




find *deduplicated.bam \
| xargs basename -s .bam \
| xargs -I{} \
${samtools} sort \
--threads ${threads} \
{}.bam \
-o {}.sorted.bam

# Index sorted files for IGV
# The "-@ ${threads}" below specifies number of CPU threads to use.

find *.sorted.bam \
| xargs basename -s .sorted.bam \
| xargs -I{} \
${samtools} index \
-@ ${threads} \
{}.sorted.bam
