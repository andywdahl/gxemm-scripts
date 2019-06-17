#!/bin/bash                         

#$ -S /bin/bash                     
#$ -cwd                            
#$ -r n                            
#$ -j y                           
#$ -l mem_free=.8G                 
#$ -l h_rt=12:00:00  

Rscript run.R
