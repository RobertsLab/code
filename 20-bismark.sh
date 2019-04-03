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


### USER NEEDS TO SET VARIABLES FOR THE FOLLOWING:
# Full path to directory with sequencing reads
reads_dir=""

# Full path to bisulftie-converted genome directory
genome_dir=""

####################################################
# DO NOT EDIT BELOW THIS LINE


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

# Concatenate R1 reads and generate lists of FastQs
for fastq in ${reads_dir}/*R1*.gz
do
  cat ${fastq} >> ${R1}
done

# Concatenate R2 reads and generate lists of FastQs
for fastq in ${reads_dir}/*R2*.gz
do
  cat ${fastq} >> ${R2}
done

## Save FastQ files to arrays
R1_array=(${reads_dir}/*_R1_*.fq*)
R2_array=(${reads_dir}/*_R2_*.fq*)

# Run bismark using bisulftie-converted genome
# Generates a set of BAM files as outputs
# Records stderr to a file for easy viewing of Bismkar summary info
${bismark_dir}/bismark \
--path_to_bowtie ${bowtie2_dir} \
--genome ${genome} \
--non_directional \
-p 28 \
-1 ${R1} \
-2 ${R2} \
2> bismark_summary.txt


# Methylation extraction
# Extracts methylation info from BAM files produced by Bismark
# Options to created a bedgraph file, counts, remove spaces from names
# and to use the "scaffolds" setting.
${bismark_dir}/bismark_methylation_extractor \
--bedgraph \
--counts \
--scaffolds \
--remove_spaces \
--multicore ${threads} \
--buffer_size 75% \
*.bam



${bismark_dir}/bismark_methylation_extractor \
--bedGraph --counts --scaffolds \
--multicore ${threads} \
--buffer_size 75% \
*deduplicated.bam



# Bismark processing report

${bismark_dir}/bismark2report

#Bismark summary report

${bismark_dir}/bismark2summary



# Sort files for methylkit and IGV

find *deduplicated.bam | \
xargs basename -s .bam | \
xargs -I{} ${samtools} \
sort --threads ${threads} {}.bam \
-o {}.sorted.bam

# Index sorted files for IGV
# The "-@ ${threads}" below specifies number of CPU threads to use.

find *.sorted.bam | \
xargs basename -s .sorted.bam | \
xargs -I{} ${samtools} \
index -@ ${threads} {}.sorted.bam
