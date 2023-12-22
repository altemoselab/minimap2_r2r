#!/bin/bash
#SBATCH --job-name=bamModEncode
#SBATCH --partition=owners,sfgf
#SBATCH --time=24:00:00
#SBATCH --mem=144G
#SBATCH -c 24
#SBATCH --chdir=/scratch/users/dubocd/altemose/slurm
#SBATCH --export=all
#SBATCH --requeue


source /home/users/dubocd/.bashrc
source activate general

sup_ref_aligned=$1
mA=$2 # aligned to sup, sorted with index 
mC=$3 # algined to sup, sorted with index

for num in {1..22} X ; do 
	samtools view -b $sup_ref_aligned chr${num} \
	| python3 /oak/stanford/groups/altemose/dubocd/miscTools/r2rModEncode.py -r $mA -a - -b A -o - --replace \
	| python3 /oak/stanford/groups/altemose/dubocd/miscTools/r2rModEncode.py -r $mC -a - -b C -o - > ${sup_ref_aligned::-3}encode.${num}.bam
done 
exit