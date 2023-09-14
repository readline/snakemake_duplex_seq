from scripts.common import (
    allocated("threads")
    ignore,
    provided 
)

# Rules
rule fastqc: 
    """
    Quality-control step to assess sequencing quality of the raw data prior removing
    adapter sequences. FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        FastQC report and zip file containing data quality information
    """

    input:
        r1 = join(workpath,"{name}.R1.fastq.gz"),
        r2 = join(workpath,"{name}.R2.fastq.gz"),
    output:
        join(workpath,"rawQC","{name}.R1_fastqc.zip"),
        join(workpath,"rawQC","{name}.R2_fastqc.zip"),
    params:
        rname  = 'rawfqc',
        outdir = join(workpath,"rawQC"),
        memory = allocated("mem", "gatk_realign", cluster).lower().rstrip('g'),
    log: 
        out = snakedir + "/logs/A1.QC/{name}.o",
        err = snakedir + "/logs/A1.QC/{name}.e"
    threads: int(allocated("threads", "fastqc_raw", cluster))
    resources: 
        mem = '16g',
        extra = ' --gres=lscratch:20 '
    container: config['images']['containername']
    shell: """
    fastqc \\
        {input.r1} \\
        {input.r2} \\
        -t {threads} \\
        -o {params.outdir}
    """