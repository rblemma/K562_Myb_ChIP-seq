#!/bin/bash

# Job name:
#SBATCH --job-name=Myb_ChIP_macs2_hg19
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
#SBATCH --mail-user=user@domain
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

Myb_Rep1=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/Myb_1_R1_NonPhiX_filt
Myb_Rep2=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/Myb_2_R1_NonPhiX_filt
Myb_Rep3=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/Myb_3_R1_NonPhiX_filt
Ctrl_Rep1=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/Ctrl_1_R1_NonPhiX_filt
Ctrl_Rep2=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/Ctrl_2_R1_NonPhiX_filt
Ctrl_Rep3=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/Ctrl_3_R1_NonPhiX_filt
inp_Myb_rep1=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/inp-Myb_1_R1_NonPhiX_filt
inp_Myb_rep2=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/inp-Myb_2_R1_NonPhiX_filt
inp_Myb_rep3=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/inp-Myb_3_R1_NonPhiX_filt
inp_Ctrl_rep1=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/inp-Ctrl_1_R1_NonPhiX_filt
inp_Ctrl_rep2=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/inp-Ctrl_2_R1_NonPhiX_filt
inp_Ctrl_rep3=/work/users/username/Roza_Myb/Myb_ChIP_hg19_SnakeMake/06_bam_filt/inp-Ctrl_3_R1_NonPhiX_filt

     
macs2 callpeak -t ${Myb_Rep1}.bam ${Myb_Rep2}.bam ${Myb_Rep3}.bam -c ${inp_Myb_rep1}.bam ${inp_Myb_rep2}.bam ${inp_Myb_rep3}.bam -f 'AUTO' -g hs -m 5 50 --bw 150 --fix-bimodal --extsize 100 --call-summits -n Myb_pooledReps_MACS2 --outdir 07_macs2_call_peaks/ -B -q 0.01 --verbose 3

macs2 callpeak -t ${Ctrl_Rep1}.bam ${Ctrl_Rep2}.bam ${Ctrl_Rep3}.bam -c ${inp_Ctrl_rep1}.bam ${inp_Ctrl_rep2}.bam ${inp_Ctrl_rep3}.bam -f 'AUTO' -g hs -m 5 50 --bw 150 --fix-bimodal --extsize 100 --call-summits -n Ctrl_pooledReps_MACS2 --outdir 07_macs2_call_peaks/ -B -q 0.01 --verbose 3
       
        
