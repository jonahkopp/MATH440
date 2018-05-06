#!/bin/bash -x
#SBATCH --job-name="Final_Proj"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=16
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --account=math540s18
#SBATCH --time=00:20:00
#SBATCH -o con.out
#SBATCH -e con.err

srun ./main_exe