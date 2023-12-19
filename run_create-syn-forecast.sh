
SLURM="#!/bin/bash\n\

#SBATCH -t 100:00:00\n\
#SBATCH --job-name=syngen\n\
#SBATCH -p normal\n\
#SBATCH --export=ALL\n\
#SBATCH --nodes=1\n\
#SBATCH --exclusive\n\
#SBATCH --mem-per-cpu=2G\n\
#SBATCH --output=syngen.txt\n\
#SBATCH --ntasks-per-node=80\n\

module load R\n\
module load gsl\n\

Rscript ./src/create_synthetic_forecasts.R"

echo -e $SLURM | sbatch 
sleep 0.5


