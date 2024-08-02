cd ~/SEM-BMI/Scripts/Shell

sex_ids=(0 1)
instances=(0 1 2 3)

for sex_id in "${sex_ids[@]}"; do
    for instance in "${instances[@]}"; do
        echo "Sex ID: $sex_id, Instance: $instance"
        
        sed "s/{sex_id}/${sex_id}/g; s/{instance}/${instance}/g" template_step2.sh > run.sh

        if [ "$sex_id" -eq 0 ]; then
            if [ "$instance" -eq 0 ]; then
                sed "/#/ s/{ncol}/14/g" run.sh > run00.sh
                qsub run00.sh
            elif [ "$instance" -eq 1 ]; then
                sed "/#/ s/{ncol}/8/g" run.sh > run01.sh
                qsub run01.sh
            elif [ "$instance" -eq 2 ]; then
                sed "/#/ s/{ncol}/11/g" run.sh > run02.sh
                qsub run02.sh
            elif [ "$instance" -eq 3 ]; then
                sed "/#/ s/{ncol}/6/g" run.sh > run03.sh
                qsub run03.sh
            else
                echo "Error occurred"
            fi
        elif [ "$sex_id" -eq 1 ]; then
            if [ "$instance" -eq 0 ]; then
                sed "/#/ s/{ncol}/12/g" run.sh > run10.sh
                qsub run10.sh
            elif [ "$instance" -eq 1 ]; then
                sed "/#/ s/{ncol}/11/g" run.sh > run11.sh
                qsub run11.sh
            elif [ "$instance" -eq 2 ]; then
                sed "/#/ s/{ncol}/8/g" run.sh > run12.sh
                qsub run12.sh
            elif [ "$instance" -eq 3 ]; then
                sed "/#/ s/{ncol}/5/g" run.sh > run13.sh
                qsub run13.sh
            else
                echo "Error occurred"
            fi
        else
            echo "Error occurred"
        fi
    done
done
