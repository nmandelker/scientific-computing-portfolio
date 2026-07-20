#!/bin/sh
#SBATCH -J Mh_Dh_15
#SBATCH -o ./logfile_%J.o
#SBATCH -e ./logfile_%J.e
#SBATCH -n 128
####SBATCH --ntasks-per-node=2
####SBATCH --ntasks=100
#SBATCH -t 3000
###SBATCH -w "exw1e2s[1-4]"
#SBATCH --partition=Apollo
#SBATCH --mail-type=end 
#SBATCH --mail-user=nir.mandelker@mail.huji.ac.il

DATE=$( date +%Y-%m-%d_%Hh%M )

##ldd ~/hello
hostname
whoami
pwd
ulimit -a


cp kh.nml 'kh'$DATE'.nml'
cd ./output
/usr/bin/mpirun --loadbalance --display-allocation \
-mca btl_openib_want_fork_support 1 \
-mca btl_openib_warn_no_device_params_found 0 \
-mca btl openib,self,sm \
-mca mpi_warn_on_fork 0 \
../ramses2d ../kh.nml > 'run_'$DATE'.log'






