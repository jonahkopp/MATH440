#!/bin/bash -x
#SBATCH --job-name="Final_Proj"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks=32
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --account=math540s18
#SBATCH --time=00:02:00
#SBATCH -o con.out
#SBATCH -e con.err

srun ./main_exe