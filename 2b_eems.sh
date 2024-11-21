#Batis

#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate eems

module load gnu/5.4.0

/common/cooperlab/phuss58/albertine/eems/runeems_snps/src/runeems_snps --params batis_aoudad_eems.params --seed 123


#Cossypha

#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate eems

module load gnu/5.4.0

/common/cooperlab/phuss58/albertine/eems/runeems_snps/src/runeems_snps --params cossypha_aoudad_eems.params --seed 123


#Chamaetylas

#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate eems

module load gnu/5.4.0

/common/cooperlab/phuss58/albertine/eems/runeems_snps/src/runeems_snps --params chamaetylas_aoudad_eems.params --seed 123


#Phylloscopus

#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate eems

module load gnu/5.4.0

/common/cooperlab/phuss58/albertine/eems/runeems_snps/src/runeems_snps --params phylloscopus_aoudad_eems.params --seed 123