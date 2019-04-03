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

#
find ${reads_dir}*_1.fq.gz \
| xargs basename -s _s1_R1_val_1.fq.gz | xargs -I{} ${bismark_dir}/bismark \
--path_to_bowtie ${bowtie2_dir} \
-genome /gscratch/srlab/sr320/data/Cvirg-genome \
-p ${threads} \
--non_directional \
-1 /gscratch/srlab/sr320/data/oakl/{}_s1_R1_val_1.fq.gz \
-2 /gscratch/srlab/sr320/data/oakl/{}_s1_R2_val_2.fq.gz \



find *.bam | \
xargs basename -s .bam | \
xargs -I{} ${bismark_dir}/deduplicate_bismark \
--bam \
--paired \
{}.bam



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
