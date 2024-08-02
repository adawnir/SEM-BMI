#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=50gb
#PBS -N clean

cd ~/SEM-BMI/Scripts

module load anaconda3/personal
source activate Renv

baseline_path=~/UKB_69328/ukb673609/Baseline/
withdrawn=/rds/general/user/rw1317/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/withdraw69328_330_20240527.txt

ukb_path=~/../projects/chadeau_ukbb_folder/live/data/project_data/UKB_673609/ukb673609.csv
repeat_path=~/UKB_69328/ukb673609/Repeated/
townsend_path=~/UKB_69328/ukb47946/extraction_townsend/outputs/ukb_extracted.rds
score_path=~/UKB_69328/Composite_scores/Data/

Rscript 0-data_preparation.R $baseline_path $withdrawn $ukb_path $repeat_path $townsend_path $score_path

