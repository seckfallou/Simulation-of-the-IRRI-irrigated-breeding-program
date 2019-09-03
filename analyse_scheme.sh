#$ -S /bin/bash
#$ -q bigmem.q
#$ -l mem_free=100G
#$ -M fallou.seck@supagro.fr
#$ -V
#$ -cwd
#$ -t 1-25

# LOAD MODULE ENVIRONMENT
module load bioinfo/R/3.5.2

Rscript --vanilla script_scheme.R

