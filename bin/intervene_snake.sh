#!/bin/bash

# Job name:
#SBATCH --job-name=Intrevene_snake
# Project:
#SBATCH --account=nn9374k
# Wall clock limit:
#SBATCH --time=4:00:00
# Max memory usage:
#SBATCH --mem-per-cpu=2000MB
# Number of cores:
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=r.b.lemma@ibv.uio.no
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

#If you need hugemem...
#--partition=hugemem

## Set up job environment:
module purge
module load python2
module load bedtools
intervene venn -i 07_macs2_call_peaks_Reps_q0.5/Myb_Rep1_MACS2_peaks.narrowPeak 07_macs2_call_peaks_Reps_q0.5/Myb_Rep2_MACS2_peaks.narrowPeak 07_macs2_call_peaks_Reps_q0.5/Myb_Rep3_MACS2_peaks.narrowPeak --names MybRep1,MybRep2,MybRep3 --bedtools-options f=0.5,r --project Myb_Rep123_intersect-three-way-venn_hg38 --save-overlaps --dpi 600 --figtype png --figsize 10 10 --fontsize 15
intervene venn -i 07_macs2_call_peaks_Reps_q0.5/Ctrl_Rep1_MACS2_peaks.narrowPeak 07_macs2_call_peaks_Reps_q0.5/Ctrl_Rep2_MACS2_peaks.narrowPeak 07_macs2_call_peaks_Reps_q0.5/Ctrl_Rep3_MACS2_peaks.narrowPeak --names CtrlRep1,CtrlRep2,CtrlRep3 --bedtools-options f=0.5,r --project Ctrl_Rep123_intersect-three-way-venn_hg38 --save-overlaps --dpi 600 --figtype png --figsize 10 10 --fontsize 15

