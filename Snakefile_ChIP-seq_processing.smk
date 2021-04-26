READ = ['Ctrl_1_R1', 'Ctrl_2_R1', 'Ctrl_3_R1', 'Myb_1_R1', 'Myb_2_R1', 'Myb_3_R1', 'inp-Ctrl_1_R1', 'inp-Ctrl_2_R1', 'inp-Ctrl_3_R1', 'inp-Myb_1_R1', 'inp-Myb_2_R1', 'inp-Myb_3_R1']
DATA='../data'


rule all:
    input:
        expand("/work/users/rozabl/Roza_Myb/Myb_ChIP_hg19_SnakeMake/01_trim-Galore/{read}.fastq.gz_trimming_report.txt", read=READ),
        expand("/work/users/rozabl/Roza_Myb/Myb_ChIP_hg19_SnakeMake/01_trim-Galore/{read}_trimmed.fq", read=READ),
        #expand("/work/users/rozabl/Roza_Myb/01_trim-Galore/{read}_trimmed.fq_fastqc.zip", read=READ),
        #expand("/work/users/rozabl/Roza_Myb/01_trim-Galore/{read}_trimmed.fq_fastqc", read=READ),
        expand("02_bbmap/{read}_NonPhiX.fq", read=READ),
        expand("03_bwa/{read}_NonPhiX.bam", read=READ),
        expand("04_bam_sort/{read}_NonPhiX_sorted.bam", read=READ),
        expand("05_bam_pcr/{read}_NonPhiX_pcr.bam", read=READ),
        expand("06_bam_filt/{read}_NonPhiX_filt.bam", read=READ),
        expand("06_bam_filt/{read}_NonPhiX_filt.bam.bai", read=READ),
        directory("07_macs2_call_peaks/"),
        expand("07_macs2_call_peak_Reps/{treatment}_{rep}_vs_{inputtreatment}_{rep}_MACS2_peaks.narrowPeak", treatment=['Myb'], inputtreatment=['inp-Myb'], rep=[1, 2, 3]),
        expand("07_macs2_call_peak_ctrl_Reps/{control}_{rep}_vs_{inputcontrol}_{rep}_MACS2_peaks.narrowPeak", control=['Ctrl'], inputcontrol=['inp-Ctrl'], rep=[1, 2, 3]),
        "08_macs2_bdgcmp/Myb_pooledReps_MACS2_bdgcmp_ppois.bdg",
        "08_macs2_bdgcmp/Ctrl_pooledReps_MACS2_bdgcmp_ppois.bdg",
        expand("09_sorted_bdg/{pkbdg}_pooledReps_MACS2_bdgcmp_ppois_sorted.bdg", pkbdg=['Myb', 'Ctrl']),
        expand("10_sorted_bw/{pkbdg}_pooledReps_MACS2_bdgcmp_ppois_sorted.bw", pkbdg=['Myb', 'Ctrl']),
        expand("08_macs2_Reps_bdgcmp/{pkbdg}_{rep}_MACS2_bdgcmp_ppois.bdg", pkbdg=['Myb', 'Ctrl'], rep=[1, 2, 3]),
        expand("10_sorted_bw/{pkbdg}_pooledReps_MACS2_bdgcmp_ppois_sorted.bw", pkbdg=['Myb', 'Ctrl']),
        expand("10_sorted_Reps_bw/{pkbdg}_{rep}_MACS2_bdgcmp_ppois_sorted.bw", pkbdg=['Myb', 'Ctrl'], rep=[1, 2, 3])
        directory("12_Intervene"),
        directory("Motif_111_MybRep1_MybRep2_MybRep3")



rule trim_galore:
    input:
        "/work/users/rozabl/Roza_Myb/Myb_ChIP_hg19_SnakeMake/{read}.fastq.gz"
    output:
        "/work/users/rozabl/Roza_Myb/Myb_ChIP_hg19_SnakeMake/01_trim-Galore/{read}.fastq.gz_trimming_report.txt",
        "/work/users/rozabl/Roza_Myb/Myb_ChIP_hg19_SnakeMake/01_trim-Galore/{read}_trimmed.fq"
        #"/work/users/rozabl/Roza_Myb/01_trim-Galore/{read}_trimmed.fq_fastqc.zip"
    params:
        out= '/work/users/rozabl/Roza_Myb/Myb_ChIP_hg19_SnakeMake/01_trim-Galore/'
    shell:
        """
        set +e
        trim_galore {input} -q 20 --phred33 --length 20 --dont_gzip --output_dir {params.out}
        set -e
        """

rule Myb_bbmap:
    input:
        "/work/users/rozabl/Roza_Myb/Myb_ChIP_hg19_SnakeMake/01_trim-Galore/{read}_trimmed.fq"
    output:
        "02_bbmap/{read}_NonPhiX.fq"
    params:
        bbmap_path= '/cluster/home/rozabl/nobackup/bbmap'
    shell:
        """
        {params.bbmap_path}/bbmap.sh ref=PhiX.fa in={input} outu={output};
        """

rule bwa_map:
    input:
        ref = "/cluster/home/rozabl/nobackup/hg19_ucsc_2bittofa/hg19.fa",
        fastq="02_bbmap/{read}_NonPhiX.fq"
    output:
        "03_bwa/{read}_NonPhiX.bam"
    threads: 15
    log:
        "logs/bwa_mem/{read}.log"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.fastq} | samtools view -bS > {output};
        """

rule samtools_sort:
    input:
        "03_bwa/{read}_NonPhiX.bam"
    output:
        "04_bam_sort/{read}_NonPhiX_sorted.bam"
    threads: 15
    shell:
        """
        samtools view -h {input} | samtools sort -n -@ {threads} | \
        samtools fixmate -r - - | samtools sort -@ {threads} | \
        samtools view -bS > {output};
        """

rule samtools_filt:
    input:
        "04_bam_sort/{read}_NonPhiX_sorted.bam"
    output:
        "05_bam_pcr/{read}_NonPhiX_pcr.bam"
    params:
        q=20
    shell:
        """
        samtools view -q {params.q} -h {input} | \
        samtools rmdup -S - - | samtools view -bS > {output};
        """

rule bam_filt_chrM:
    input:
        "05_bam_pcr/{read}_NonPhiX_pcr.bam"
    output:
        "06_bam_filt/{read}_NonPhiX_filt.bam"
    shell:
        """
        samtools view -h {input} | awk '{{if($3!="chrM")print}}' | \
        samtools view -bS > {output};
        """

rule samtools_index:
    input:
        "06_bam_filt/{read}_NonPhiX_filt.bam"
    output:
        "06_bam_filt/{read}_NonPhiX_filt.bam.bai"
    shell:
        """
        samtools index {input} > {output};
        """

rule macs2_call_peaks_pooledReps:
    input:
        trtreps=expand("06_bam_filt/{treatment}_{rep}_R1_NonPhiX_filt.bam", treatment=['Myb'], rep=[1, 2, 3]),
        ctrlreps=expand("06_bam_filt/{control}_{rep}_R1_NonPhiX_filt.bam", control=['Ctrl'], rep=[1, 2, 3]),
        intrtreps=expand("06_bam_filt/{inputtreatment}_{rep}_R1_NonPhiX_filt.bam", inputtreatment=['inp-Myb'], rep=[1, 2, 3]),
        inctrlreps=expand("06_bam_filt/{inputcontrol}_{rep}_R1_NonPhiX_filt.bam", inputcontrol=['inp-Ctrl'], rep=[1, 2, 3])

    output:
        directory("07_macs2_call_peaks/")

    params:
        bw=150,
        es=100,
        od="07_macs2_call_peaks/",
        qv=0.01,
        vb=3,
        file='AUTO'

    shell:
        """
        conda create --name py2 python=2.7;
        conda activate py2;
        module load macs2;
        macs2 callpeak -t {input.trtreps} -c {input.intrtreps} -f {params.file} -g hs \
            -m 5 50 --bw {params.bw} --fix-bimodal --extsize {params.es} --call-summits \
            -n Myb_pooledReps_MACS2 --outdir {params.od} -B -q {params.qv} --verbose {params.vb};
        macs2 callpeak -t {input.ctrlreps} -c {input.inctrlreps} -f {params.file} -g hs \
            -m 5 50 --bw {params.bw} --fix-bimodal --extsize {params.es} --call-summits \
            -n Ctrl_pooledReps_MACS2 --outdir {params.od} -B -q {params.qv} --verbose {params.vb};
        module unload macs2;
        conda deactivate;
        """

rule macs2_bdgcmp_refine_peaks:
    input:
        trt_pileup= "07_macs2_call_peaks/Myb_pooledReps_MACS2_treat_pileup.bdg",
        trt_lambda= "07_macs2_call_peaks/Myb_pooledReps_MACS2_control_lambda.bdg",
        ctrl_pileup= "07_macs2_call_peak_ctrl_Reps/Ctrl_pooledReps_MACS2_treat_pileup.bdg",
        ctrl_lambda= "07_macs2_call_peak_ctrl_Reps/Ctrl_pooledReps_MACS2_control_lambda.bdg"

    output:
        "08_macs2_bdgcmp/Myb_pooledReps_MACS2_bdgcmp_ppois.bdg",
        "08_macs2_bdgcmp/Ctrl_pooledReps_MACS2_bdgcmp_ppois.bdg"
    params:
        od="08_macs2_bdgcmp/"

    shell:
        """
        conda activate py2;
        module load macs2;
        macs2 bdgcmp -t {input.trt_pileup} -c {input.trt_lambda} -m ppois \
            --outdir {params.od} --o-prefix Myb_pooledReps_MACS2_bdgcmp;
        macs2 bdgcmp -t {input.ctrl_pileup} -c {input.ctrl_lambda} -m ppois \
            --outdir {params.od} --o-prefix Ctrl_pooledReps_MACS2_bdgcmp;
        module unload macs2;
        conda deactivate;
        """

rule sort_bedgraph:
    input:
        expand("08_macs2_bdgcmp/{{pkbdg}}_pooledReps_MACS2_bdgcmp_ppois.bdg", pkbdg=['Myb', 'Ctrl'])
    output:
        "09_sorted_bdg/{pkbdg}_pooledReps_MACS2_bdgcmp_ppois_sorted.bdg"
    shell:
        """
        sort -k1,1 -k2,2n {input} > {output};
        """

rule bedgraph2bigwig:
    input:
        bedgraphs=expand("09_sorted_bdg/{{pkbdg}}_pooledReps_MACS2_bdgcmp_ppois_sorted.bdg", pkbdg=['Myb', 'Ctrl']),
        chromsizes="{0}/UCSC/hg19.chrom.sizes".format(DATA)
    output:
        "10_sorted_bw/{pkbdg}_pooledReps_MACS2_bdgcmp_ppois_sorted.bw"
    shell:
        """
        ../bin/bedGraphToBigWig {input.bedgraphs} {input.chromsizes} {output};
        """

#### For individual Reps ##
rule macs2_call_peaks_MyB_Reps:
    input:
        trtreps=expand("06_bam_filt/{{treatment}}_{{rep}}_R1_NonPhiX_filt.bam", treatment=['Myb'], rep=[1, 2, 3]),
        intrtreps=expand("06_bam_filt/{{inputtreatment}}_{{rep}}_R1_NonPhiX_filt.bam", inputtreatment=['inp-Myb'], rep=[1, 2, 3]),

    output:
        "07_macs2_call_peak_Reps/{treatment}_{rep}_vs_{inputtreatment}_{rep}_MACS2_peaks.narrowPeak",

    params:
        bw=150,
        es=100,
        od="07_macs2_call_peak_Reps/",
        qv=0.01,
        vb=3,
        file='AUTO',
        no_myb=expand('{{treatment}}_{{rep}}_vs_{{inputtreatment}}_{{rep}}', treatment=['Myb'], inputtreatment=['inp-Myb'], rep=[1, 2, 3])

    shell:
        """
        conda create --name py2 python=2.7;
        conda activate py2;
        module load macs2;
        macs2 callpeak -t {input.trtreps} -c {input.intrtreps} -f {params.file} -g hs \
            -m 5 50 --bw {params.bw} --fix-bimodal --extsize {params.es} --call-summits \
            -n {params.no_myb}_MACS2 --outdir {params.od} -B -q {params.qv} --verbose {params.vb};
        conda deactivate;
        """

rule macs2_call_peaks_Ctrl_Reps:
    input:
        ctrlreps=expand("06_bam_filt/{{control}}_{{rep}}_R1_NonPhiX_filt.bam", control=['Ctrl'], rep=[1, 2, 3]),
        inctrlreps=expand("06_bam_filt/{{inputcontrol}}_{{rep}}_R1_NonPhiX_filt.bam", inputcontrol=['inp-Ctrl'], rep=[1, 2, 3])

    output:
        "07_macs2_call_peak_ctrl_Reps/{control}_{rep}_vs_{inputcontrol}_{rep}_MACS2_peaks.narrowPeak"

    params:
        bw=150,
        es=100,
        od="07_macs2_call_peak_ctrl_Reps/",
        qv=0.01,
        vb=3,
        file='AUTO',
        no_ctrl=expand('{{control}}_{{rep}}_vs_{{inputcontrol}}_{{rep}}', treatment=['Ctrl'], inputtreatment=['inp-Ctrl'], rep=[1, 2, 3])

    shell:
        """
        conda create --name py2 python=2.7;
        conda activate py2;
        module load macs2;
        macs2 callpeak -t {input.ctrlreps} -c {input.inctrlreps} -f {params.file} -g hs \
            -m 5 50 --bw {params.bw} --fix-bimodal --extsize {params.es} --call-summits \
            -n {params.no_ctrl}_MACS2 --outdir {params.od} -B -q {params.qv} --verbose {params.vb};
        conda deactivate;
        """

rule macs2_bdgcmp_refine_peaks_Reps:
    input:
        trt_pileup= expand("07_macs2_call_peak_Reps/{treatment}_{rep}_vs_{inputtreatment}_{rep}_MACS2_treat_pileup.bdg",
            treatment=['Myb'], inputtreatment=['inp-Myb'], rep=[1, 2, 3]),
        trt_lambda= expand("07_macs2_call_peak_Reps/{treatment}_{rep}_vs_{inputtreatment}_{rep}_MACS2_control_lambda.bdg",
            treatment=['Myb'], inputtreatment=['inp-Myb'], rep=[1, 2, 3]),
        ctrl_pileup= expand("07_macs2_call_peak_ctrl_Reps/{control}_{rep}_vs_{inputcontrol}_{rep}_MACS2_treat_pileup.bdg",
            control=['Ctrl'], inputcontrol=['inp-Ctrl'], rep=[1, 2, 3]),
        ctrl_lambda= expand("07_macs2_call_peak_ctrl_Reps/{control}_{rep}_vs_{inputcontrol}_{rep}_MACS2_control_lambda.bdg",
            control=['Ctrl'], inputcontrol=['inp-Ctrl'], rep=[1, 2, 3])

    output:
        #"08_macs2_Reps_bdgcmp/{pkbdg}_{rep}_MACS2_bdgcmp_ppois.bdg",
        expand("08_macs2_Reps_bdgcmp/{pkbdg}_{rep}_MACS2_bdgcmp_ppois.bdg", pkbdg=['Myb', 'Ctrl'], rep=[1, 2, 3])
    params:
        od="08_macs2_bdgcmp/",
        pref_myb=expand("{pkbdg}_{rep}_MACS2_bdgcmp_ppois", pkbdg=['Myb'], rep=[1, 2, 3]),
        pref_ctrl=expand("{pkbdg}_{rep}_MACS2_bdgcmp_ppois", pkbdg=['Ctrl'], rep=[1, 2, 3])

    shell:
        """
        conda activate py2;
        module load macs2;
        macs2 bdgcmp -t {input.trt_pileup} -c {input.trt_lambda} -m ppois \
            --outdir {params.od} --o-prefix {params.pref_myb};
        macs2 bdgcmp -t {input.ctrl_pileup} -c {input.ctrl_lambda} -m ppois \
            --outdir {params.od} --o-prefix {params.pref_ctrl};
        module unload macs2;
        conda deactivate;
        """

rule sort_bedgraph_Reps:
    input:
        expand("08_macs2_Reps_bdgcmp/{pkbdg}_{rep}_MACS2_bdgcmp_ppois.bdg", pkbdg=['Myb', 'Ctrl'], rep=[1, 2, 3])
    output:
        "09_sorted_Reps_bdg/{pkbdg}_{rep}_MACS2_bdgcmp_ppois_sorted.bdg"
    shell:
        """
        sort -k1,1 -k2,2n {input} > {output};
        """

rule bedgraph2bigwig_Reps:
    input:
        bedgraphs= expand("09_sorted_Reps_bdg/{pkbdg}_{rep}_MACS2_bdgcmp_ppois_sorted.bdg", pkbdg=['Myb', 'Ctrl'], rep=[1, 2, 3]),
        chromsizes="{0}/UCSC/hg19.chrom.sizes".format(DATA)
    output:
        "10_sorted_Reps_bw/{pkbdg}_{rep}_MACS2_bdgcmp_ppois_sorted.bw"
    shell:
        """
        ../bin/bedGraphToBigWig {input.bedgraphs} {input.chromsizes} {output};
        """


rule bedtools_merge:
    input:
        myb_beds=expand("07_macs2_call_peak_Reps/{{treatment}}_{{rep}}_vs_{{inputtreatment}}_{{rep}}_MACS2_peaks.narrowPeak",
            treatment=['Myb'], inputtreatment=['inp-Myb'], rep=[1, 2, 3]),
        ctrl_beds=expand("07_macs2_call_peak_ctrl_Reps/{{control}}_{{rep}}_vs_{{inputcontrol}}_{{rep}}_MACS2_peaks.narrowPeak",
            control=['Ctrl'], inputcontrol=['inp-Ctrl'], rep=[1, 2, 3])
    output:
        myb_beds=expand("07_macs2_call_peak_Reps/{treatment}_{rep}_vs_{inputtreatment}_{rep}_MACS2_peaks.narrowPeak_merged.bed",
            treatment=['Myb'], inputtreatment=['inp-Myb'], rep=[1, 2, 3]),
        ctrl_beds=expand("07_macs2_call_peak_ctrl_Reps/{control}_{rep}_vs_{inputcontrol}_{rep}_MACS2_peaks.narrowPeak_merged.bed",
            control=['Ctrl'], inputcontrol=['inp-Ctrl'], rep=[1, 2, 3])
    shell:
        """
        modulel oad bedtools;
        bedtools merge -i <(sort -k1,1 -k2,2n {input.myb_beds}) > {output.myb_beds};
        bedtools merge -i <(sort -k1,1 -k2,2n {input.ctrl_beds}) > {output.ctrl_beds};
        module unload bedtools;
        """

rule intervene:
    input:
        myb_beds="07_macs2_call_peak_Reps/{treatment}_{rep}_vs_{inputtreatment}_{rep}_MACS2_peaks.narrowPeak_merged.bed",
        ctrl_beds="07_macs2_call_peak_ctrl_Reps/{control}_{rep}_vs_{inputcontrol}_{rep}_MACS2_peaks.narrowPeak_merged.bed"
    output:
        directory("12_Intervene")
    shell:
        """
        module load bedtools;
        intervene venn -i {input.myb_beds} --names MybRep1,MybRep2,MybRep3 --bedtools-options f=0.5,r \
            --project Myb_Rep123_bedtools_merge_intersect --save-overlaps --dpi 600 --figtype png \
            --figsize 10 10 --fontsize 15 --colors darkturquoise,deepskyblue,mediumpurple;
        intervene venn -i {input.ctrl_beds} --names CtrlRep1,CtrlRep2,CtrlRep3 --bedtools-options f=0.5,r \
            --project Ctrl_Rep123_bedtools_merge_intersect --save-overlaps --dpi 600 --figtype png \
            --figsize 10 10 --fontsize 15 --colors darkturquoise,deepskyblue,mediumpurple;
        """

rule HOMER_motif:
    input:
        "12_Intervene/Intervene_results/sets/111_MybRep1_MybRep2_MybRep3.ed"
    output:
        directory("Motif_111_MybRep1_MybRep2_MybRep3")
    shell:
        """
        findMotifsGenome.pl {input} hg19 {output} -size given -preparsedDir 12_Intervene/Intervene_results/sets/;
        """
