#!/bin/bash

#SBATCH --account=psych
#SBATCH --job-name=1st_level_conf
#SBATCH -c 24
#SBATCH --time=11:55:00
#SBATCH --mem=120gb

bash /rigel/psych/users/cmh2228/dynCon/rl_flexibility/1st_level_conf.sh $arg1 $arg2 $arg3 $arg4 $arg5
