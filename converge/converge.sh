#!/bin/bash                         

#$ -S /bin/bash                     
#$ -cwd                            
#$ -r n                            
#$ -j y                           
#$ -l mem_free=8G
#$ -l h_rt=42:00:00  

Rscript run.R
