#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=8:mem=30gb
#PBS -N step1
#PBS -J 0-1
#PBS -o /rds/general/user/rw1317/home/SEM-BMI/Jobs/
#PBS -e /rds/general/user/rw1317/home/SEM-BMI/Jobs/

cd ~/SEM-BMI/Scripts

module load anaconda3/personal
source activate Renv

Rscript 1-direct_effects.R $PBS_ARRAY_INDEX {instance}



