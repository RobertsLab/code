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



# Directories and programs
wd=$(pwd)
bismark_dir="/gscratch/srlab/programs/Bismark-0.21.0"
bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/"
samtools="/gscratch/srlab/programs/samtools-1.9/samtools"
reads_dir="/gscratch/srlab/sr320/data/oakl/"


# what is this line for ?
source /gscratch/srlab/programs/scripts/paths.sh



#
find ${reads_dir}*_1.fq.gz \
| xargs basename -s _s1_R1_val_1.fq.gz | xargs -I{} ${bismark_dir}/bismark \
--path_to_bowtie ${bowtie2_dir} \
-genome /gscratch/srlab/sr320/data/Cvirg-genome \
-p 14 \
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
--multicore 14 \
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
sort --threads 28 {}.bam \
-o {}.sorted.bam

# Index sorted files for IGV
# The "-@ 16" below specifies number of CPU threads to use.

find *.sorted.bam | \
xargs basename -s .sorted.bam | \
xargs -I{} ${samtools} \
index -@ 28 {}.sorted.bam
