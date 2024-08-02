cd ~/SEM-BMI/Scripts/Shell

for i in $(seq 0 1 3); do
echo $i

sed "s/{instance}/${i}/g" template_step1.sh > run.sh

qsub run.sh

done