#!/bin/bash
#SBATCH --partition=herringlab-low,statdept-low
#SBATCH --account=statdept
#SBATCH -c1

singularity run -B /work/sta790/ff31:/work/sta790/ff31 felpo_ex.sif factor_interactions/codes/cluster_longeneker.R 4000 3000 1

