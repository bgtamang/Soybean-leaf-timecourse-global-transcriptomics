#!/bin/bash 
#SBATCH -N 1 
#SBATCH -n 4
#SBATCH --mem=64G 
#SBATCH --mail-user=bgtamang@illinois.edu 
#SBATCH --mail-type=ALL
#SBATCH --job-name=salmon_map_60 
#SBATCH -D /home/a-m/bgtamang/soybean_RNASeq_2025/src/slurm-out 
#SBATCH -p normal
#SBATCH --array=1-60%4

# NOTES on options above: 1. If you are re-using this script, please change the --mail-user above to your e-mail 
# address!!! 2. The -D location was set primarily to control where the slurm-*.out files end up
#    FULL PATH NAME MUST BE USED HERE!! 3. The -p normal partition is required for normal alignment - cannot be done 
# on classroom partition if have
#    normal-sized (>= 10M) fastq files 4. On normal sized files and genomes, you will need more threads and more 
# memory. change the directory to the main one containing data/, results/ and src/. All relative file paths in the 
# codes below will be to this location

cd /home/a-m/bgtamang/soybean_RNASeq_2025/

#load the program needed, name may change in future so check each time with 'module avail' on command line

module load Salmon/1.10.0-IGB-gcc-8.2.0

#set up variable to pull base file names - get this from our STAR or featureCounts scripts??
line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" src/basenamesFULL.txt)

#make directory if it does not exist

mkdir -p results/salmon

#run salmon

echo "start Salmon mapping"
     salmon quant -i data/genome/Salmon-1.10.0-Soybean_Index -l A \
    -r data/raw-seq/${line}.fastq.gz \
    --numBootstraps=30 \
    -o results/salmon/${line} -p $SLURM_NTASKS \
    --validateMappings \
    --gcBias \
    --seqBias 
echo "end Salmon mapping"
