#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o Rout/out			#-- output directory (fill in)
#$ -e Rout/err                     #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r n                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=7G                  #-- submits on nodes with enough free memory (required)
#$ -l h_rt=8:20:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1-200													#-- remove first '#' to specify the number of

date
hostname

Rscript run.R $SGE_TASK_ID
