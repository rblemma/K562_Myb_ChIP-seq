#!/bin/bash

# Job name:
#SBATCH --job-name=DU145_ATAC_bwa
#
# Project:
#SBATCH --account=nn9632k
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=6000M
#
#SBATCH --cpus-per-task=10

# notify job failure
#SBATCH --mail-user=rozabl@ncmm.uio.no
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

#If you need hugemem...
#--partition=hugemem

## Set up the job environment
source /cluster/bin/jobsetup
set -o errexit

##Load job
module purge
module load macs2


macs2 bdgcmp -t 07_macs2_call_peaks/Myb_pooledReps_MACS2_treat_pileup.bdg -c 07_macs2_call_peaks/Myb_pooledReps_MACS2_control_lambda.bdg -m ppois --outdir ./08_macs2_bdgcmp/ --o-prefix Myb_pooledReps_MACS2_bdgcmp

macs2 bdgcmp -t 07_macs2_call_peaks/Ctrl_pooledReps_MACS2_treat_pileup.bdg -c 07_macs2_call_peaks/Ctrl_pooledReps_MACS2_control_lambda.bdg -m ppois --outdir ./08_macs2_bdgcmp/ --o-prefix Ctrl_pooledReps_MACS2_bdgcmp
       
        
