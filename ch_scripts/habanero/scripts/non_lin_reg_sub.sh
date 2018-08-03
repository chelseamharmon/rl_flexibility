#!/bin/bash

#SBATCH --account=psych
#SBATCH --job-name=non_lin_reg
#SBATCH -c 24
#SBATCH --time=11:55:00
#SBATCH --mem=120gb

bash /rigel/psych/users/cmh2228/dynCon/scripts/non_lin_reg.sh $arg1 $arg2


