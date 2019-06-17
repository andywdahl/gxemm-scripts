#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -o Rout_HE/mainout			#-- output directory (fill in)
#$ -e Rout_HE/err                     #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -l mem_free=.5G                  #-- submits on nodes with enough free memory (required)
#$ -l h_rt=12:20:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1:500													#-- remove first '#' to specify the number of

date
hostname

Rscript run_HE.R $SGE_TASK_ID
