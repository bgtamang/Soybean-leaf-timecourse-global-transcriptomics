#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bgtamang@illinois.edu
#SBATCH -J salmon.prep
#SBATCH -D /home/a-m/bgtamang/soybean_RNASeq_2025/src/slurm-out
#SBATCH -p normal

# See https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/


# Change to the directory with the reference genome

cd ~/soybean_RNASeq_2025/data/genome

# Make a text file listing the names of the "decoys".

grep "^>" <(gunzip -c Gmax_880_v6.0.fa.gz) | cut -d " " -f 1 > decoys.txt 
sed -i.bak -e 's/>//g' decoys.txt

# Combine the transcriptome and genome sequences into one file
cat Gmax_880_Wm82.a6.v1.transcript.fa.gz Gmax_880_v6.0.fa.gz > gentrome.fa.gz
