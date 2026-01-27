#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bgtamang@illinois.edu
#SBATCH -J salmon.index
#SBATCH -D /home/a-m/bgtamang/soybean_RNASeq_2025/src/slurm-out
#SBATCH -p normal

cd /home/a-m/bgtamang/soybean_RNASeq_2025

module load Salmon/1.10.0-IGB-gcc-8.2.0

salmon index \
 -t ~/soybean_RNASeq_2025/data/genome/gentrome.fa.gz \
 -i ~/soybean_RNASeq_2025/data/genome/Salmon-1.10.0-Soybean_Index \
 -d ~/soybean_RNASeq_2025/data/genome/decoys.txt \
 -p $SLURM_NTASKS

echo "Yay..done"
