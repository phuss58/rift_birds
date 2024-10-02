#!/bin/bash
#SBATCH --err=%j.err
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-39

# load bcftools
module load vcftools

# pull from /common/
workdir=/common/cooperlab/phuss58/albertine/

# save to /work/
savedir=/work/cooperlab/phuss58/albertine/

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep ${workdir}/batis.txt --max-missing 1.0 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --mac 2 --recode --recode-INFO-all --remove-indels --out ${savedir}/05_filtered_vcf/batis_structure_${region_array}_structure_nowindow
