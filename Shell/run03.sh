#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=8:mem=30gb
#PBS -N step2
#PBS -J 1-6
#PBS -o /rds/general/user/rw1317/home/SEM-BMI/Jobs/
#PBS -e /rds/general/user/rw1317/home/SEM-BMI/Jobs/

cd ~/SEM-BMI/Scripts

module load anaconda3/personal
source activate Renv

Rscript 2-indirect_effects.R $PBS_ARRAY_INDEX 0 3

