#!/bin/bash
#SBATCH -p skx-normal
#SBATCH -A TG-AST190003
#SBATCH --job-name=PT_134_151
#SBATCH -N 1
#SBATCH -n 32
##SBATCH --mem=30000
#SBATCH --time=0:25:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nir_mandelker@ucsb.edu
set -e

module purge
module load intel/18.0.2  impi/18.0.2
module load launcher/3.1

launcher_file="launcher_file_ART_loop.txt"

date

rm -rf $launcher_file
touch $launcher_file
for i in {134..151}; do
   echo "./loop_output_ART.exe $i $i" >> $launcher_file
done

export LAUNCHER_PPN=20
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=$launcher_file

$LAUNCHER_DIR/paramrun

#DATE=$( date +%Y-%m-%d_%Hh%M )

#hostname
#whoami
#pwd
#ulimit -a
