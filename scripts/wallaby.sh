#!/bin/bash
#
#SBATCH --job-name=stingray
#SBATCH --ntasks=6
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=16000

module load gfortran/6.3.0 hdf5/1.10.2

for seed in {1..5}
do
   srun -n 1 /home/dobreschkow/stingray/stingray -parameterset wallaby-micro-hyades_$seed -parameterfile /home/dobreschkow/stingray/parameters.txt -logfile /home/dobreschkow/log_wallaby_micro_$seed.txt &
   srun -n 1 /home/dobreschkow/stingray/stingray -parameterset wallaby-medi-hyades_$seed -parameterfile /home/dobreschkow/stingray/parameters.txt -logfile /home/dobreschkow/log_wallaby_medi_$seed.txt &
done
wait