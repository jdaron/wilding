#!/bin/bash 
#SBATCH --job-name=qualimap ## On définit le nom du job
#SBATCH -n 1 ## On définit le nombre de tâches
#SBATCH -c 1 ## On définit le nombre cpu
#SBATCH --export variables ## On exporte les variables d'environnement
#SBATCH -p normal
# #SBATCH --nodelist=node14

path_to_tmp="/scratch/daron_$SLURM_JOB_ID"
echo $path_to_tmp
mkdir $path_to_tmp/
path_to_dir=($PWD)
inputFile=${1}

scp nas3:/data3/projects/plasmodium/anopheles/bam/$inputFile.bam $path_to_tmp/

cd $path_to_tmp/
mkdir out

qualimap bamqc -bam $inputFile.bam -c --java-mem-size=8G -outdir out --outfile $inputFile.bamqc -nt 1 -outformat PDF

mv out/genome_results.txt out/$inputFile.genome_results.txt
scp out/* nas3:$path_to_dir/.


rm -rf $path_to_tmp
