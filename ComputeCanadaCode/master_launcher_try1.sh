#!/bin/bash
for i in `seq 70 96`;
do

sbatch --export=my_var=$i --job-name=randomrep_try1_$i --output=output_randomreptry1/output_randomrep_try1_$i.txt --error=error_randomrep_try1/error_randomreptry1_$i.txt randomreptry1.sh
done