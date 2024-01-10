import os,sys
import yaml
from os.path import join
import pandas as pd
from scripts.load import samplesheet
from scripts.utils import (allocated, ignore)

# Load config file
configfile: 'config.yaml'
workdir: config['workdir']

pipedir = config['pipelinedir']
print(config)

# Load cluster config
with open(join(config['pipelinedir'], 'cluster.yaml')) as infile:
    cluster = yaml.safe_load(infile)

# Load sample sheet
sampledic, rundic, run2sample = samplesheet(join(config['pipelinedir'],'samplesheet.tsv'))

print(sampledic)
print(rundic)
print(run2sample)


rule all:
    input: 
        expand(
            join(config['workdir'], "01.MergeData", "{sample}", "{sample}.QC.html"),
            sample=sampledic,
        ),

        expand(
            join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "{sample}.report.html"),
            sample=sampledic,
        ),

        expand(
            join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "dcs", "{sample}.dcs.snps.anno.vcf.gz"),
            sample=sampledic,
        ),

        expand(
            join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "dcs", "{sample}.dcs.anno.vcf.gz"),
            sample=sampledic,
        ),

        expand(
            join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "sscs", "{sample}.sscs.snps.anno.vcf.gz"),
            sample=sampledic,
        ),

        expand(
            join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "sscs", "{sample}.sscs.anno.vcf.gz"),
            sample=sampledic,
        ),

        report  = join(config['workdir'], "03.MultiQC", "multiqc_report.html"),
        summary = join(config['workdir'], "All.Duplex.summary.tsv"),

rule Merge1:
    input:
        reads1 = lambda wildcards: [rundic[rg]['r1'] for rg in sampledic[wildcards.sample]],
    output:
        reads1out = temp(join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R1.fq.gz")),
    log: 
        out = join(config['pipelinedir'], "logs", "Merge1", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Merge1", "{sample}.e"),
    threads:
        int(allocated("threads", "Merge1", cluster))
    container:
        config['container']['duplex']
    shell:
        """
        if [ $(ls -1 {input.reads1} | wc -l) -eq 1 ]; then
            ln -s {input.reads1} {output.reads1out} \
               >> {log.out} 2> {log.err}
        else
            zcat {input.reads1} | gzip -c - > {output.reads1out} 2> {log.err}
        fi
        """

rule Merge2:
    input:
        reads2 = lambda wildcards: [rundic[rg]['r2'] for rg in sampledic[wildcards.sample]],
    output:
        reads2out = temp(join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R2.fq.gz")),
    log: 
        out = join(config['pipelinedir'], "logs", "Merge2", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Merge2", "{sample}.e"),
    threads:
        int(allocated("threads", "Merge2", cluster))
    shell:
        """
        if [ $(ls -1 {input.reads2} | wc -l) -eq 1 ]; then
            ln -s {input.reads2} {output.reads2out} \
               >> {log.out} 2> {log.err}
        else
            zcat {input.reads2} | gzip -c - > {output.reads2out} 2> {log.err}
        fi
        """

rule Fastp:
    input:
        reads1 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R1.fq.gz"),
        reads2 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R2.fq.gz"),
    output:
        reads1 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R1.cln.fq.gz"),
        reads2 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R2.cln.fq.gz"),
        htmlout = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.QC.html"),
        jsonout = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.QC.json"),
    log:
        out = join(config['pipelinedir'], "logs", "Fastp", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Fastp", "{sample}.e"),
    threads:
        int(allocated("threads", "Fastp", cluster))
    container:
        config['container']['duplex']
    shell:
        '''
        fastp -i {input.reads1} \
            -I {input.reads2} \
            -o {output.reads1} \
            -O {output.reads2} \
            -h {output.htmlout} \
            -j {output.jsonout} \
            -w {threads} \
            >> /dev/null 2> {log.err}
        '''

rule Duplex:
    input:
        reads1 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R1.cln.fq.gz"),
        reads2 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R2.cln.fq.gz"),
    output:
        report = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "{sample}.report.html"),
        summary= join(config['workdir'], "02.Duplex", "{sample}", "{sample}.config.csv.summary.csv"),
        dcssnp = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "dcs","{sample}.dcs.snps.vcf"),
        dcssnv = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "dcs","{sample}.dcs.vcf"),
        sscssnp = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "sscs","{sample}.sscs.snps.vcf"),
        sscssnv = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "sscs","{sample}.sscs.vcf"),
    params:
        path   = join(config['workdir'], "02.Duplex","{sample}","workdir"),
        config = join(config['pipelinedir'], "duplex_config", "{sample}.config.csv"),
    log:
        out = join(config['pipelinedir'], "logs", "Duplex", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Duplex", "{sample}.e"),
    threads:
        int(allocated("threads", "Duplex", cluster))
    container:
        config['container']['duplex']
    shell:
        '''
        mkdir -p {params.path} && \
        cd {params.path} && \
        cp {params.config} ../ && \
        ln -s {input.reads1} read1.fq.gz && \
        ln -s {input.reads2} read2.fq.gz && \
        cd .. && \
        DS {wildcards.sample}.config.csv > {log.out} 2>{log.err} 
        '''

rule Annotation1:
    input:
        dcssnp = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "dcs","{sample}.dcs.snps.vcf"),
    output:
        dcssnp = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "dcs","{sample}.dcs.snps.anno.vcf.gz"),
    params:
        ref = join(config['cachedir'], "reference","snpeff"),
        dbnsfp = join(config['cachedir'], "reference","dbNSFP4.4a_hg38.txt.gz"),
        gnomad = join(config['cachedir'], "reference","gnomad"),
    log:
        out = join(config['pipelinedir'], "logs", "Annotation1", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Annotation1", "{sample}.e"),
    threads:
        int(allocated("threads", "Annotation", cluster))
    container:
        config['container']['duplex']
    shell:
        """
a={input.dcssnp}
b={output.dcssnp}
snpEff -Xmx14g -noStats -dataDir {params.ref} -lof -motif -hgvs GRCh38.105 $a > $a.temp.vcf 2>{log.err}
SnpSift -Xmx14g dbnsfp -db {params.dbnsfp} $a.temp.vcf >$a.temp.vcf2 2>>{log.err} 
mv $a.temp.vcf2 $a.temp.vcf
for chrid in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
    SnpSift -Xmx14g annotate -info AF,AF_nfe,grpmax,fafmax_faf95_max,nhomalt {params.gnomad}/gnomad.genomes.v4.0.sites.chr$chrid.vcf.bgz $a.temp.vcf > $a.temp.vcf2 2>>{log.err}
    mv $a.temp.vcf2 $a.temp.vcf
done
cat $a.temp.vcf|bgzip -c - > $b
tabix -p vcf $b
"""

rule Annotation2:
    input:
        dcssnv = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "dcs","{sample}.dcs.vcf"),
    output:
        dcssnv = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "dcs","{sample}.dcs.anno.vcf.gz"),
    params:
        ref = join(config['cachedir'], "reference","snpeff"),
        dbnsfp = join(config['cachedir'], "reference","dbNSFP4.4a_hg38.txt.gz"),
        gnomad = join(config['cachedir'], "reference","gnomad"),
    log:
        out = join(config['pipelinedir'], "logs", "Annotation2", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Annotation2", "{sample}.e"),
    threads:
        int(allocated("threads", "Annotation", cluster))
    container:
        config['container']['duplex']
    shell:
        """
a={input.dcssnv}
b={output.dcssnv}
snpEff -Xmx14g -noStats -dataDir {params.ref} -lof -motif -hgvs GRCh38.105 $a > $a.temp.vcf 2>{log.err}
SnpSift -Xmx14g dbnsfp -db {params.dbnsfp} $a.temp.vcf >$a.temp.vcf2 2>>{log.err} 
mv $a.temp.vcf2 $a.temp.vcf
for chrid in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
    SnpSift -Xmx14g annotate -info AF,AF_nfe,grpmax,fafmax_faf95_max,nhomalt {params.gnomad}/gnomad.genomes.v4.0.sites.chr$chrid.vcf.bgz $a.temp.vcf > $a.temp.vcf2 2>>{log.err}
    mv $a.temp.vcf2 $a.temp.vcf
done
cat $a.temp.vcf|bgzip -c - > $b
tabix -p vcf $b
"""

rule Annotation3:
    input:
        sscssnp = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "sscs","{sample}.sscs.snps.vcf"),
    output:
        sscssnp = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "sscs","{sample}.sscs.snps.anno.vcf.gz"),
    params:
        ref = join(config['cachedir'], "reference","snpeff"),
        dbnsfp = join(config['cachedir'], "reference","dbNSFP4.4a_hg38.txt.gz"),
        gnomad = join(config['cachedir'], "reference","gnomad"),
    log:
        out = join(config['pipelinedir'], "logs", "Annotation3", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Annotation3", "{sample}.e"),
    threads:
        int(allocated("threads", "Annotation", cluster))
    container:
        config['container']['duplex']
    shell:
        """
a={input.sscssnp}
b={output.sscssnp}
snpEff -Xmx14g -noStats -dataDir {params.ref} -lof -motif -hgvs GRCh38.105 $a > $a.temp.vcf 2>{log.err}
SnpSift -Xmx14g dbnsfp -db {params.dbnsfp} $a.temp.vcf >$a.temp.vcf2 2>>{log.err} 
mv $a.temp.vcf2 $a.temp.vcf
for chrid in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
    SnpSift -Xmx14g annotate -info AF,AF_nfe,grpmax,fafmax_faf95_max,nhomalt {params.gnomad}/gnomad.genomes.v4.0.sites.chr$chrid.vcf.bgz $a.temp.vcf > $a.temp.vcf2 2>>{log.err}
    mv $a.temp.vcf2 $a.temp.vcf
done
cat $a.temp.vcf|bgzip -c - > $b
tabix -p vcf $b
"""

rule Annotation4:
    input:
        sscssnv = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "sscs","{sample}.sscs.vcf"),
    output:
        sscssnv = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "sscs","{sample}.sscs.anno.vcf.gz"),
    params:
        ref = join(config['cachedir'], "reference","snpeff"),
        dbnsfp = join(config['cachedir'], "reference","dbNSFP4.4a_hg38.txt.gz"),
        gnomad = join(config['cachedir'], "reference","gnomad"),
    log:
        out = join(config['pipelinedir'], "logs", "Annotation4", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Annotation4", "{sample}.e"),
    threads:
        int(allocated("threads", "Annotation", cluster))
    container:
        config['container']['duplex']
    shell:
        """
a={input.sscssnv}
b={output.sscssnv}
snpEff -Xmx14g -noStats -dataDir {params.ref} -lof -motif -hgvs GRCh38.105 $a > $a.temp.vcf 2>{log.err}
SnpSift -Xmx14g dbnsfp -db {params.dbnsfp} $a.temp.vcf >$a.temp.vcf2 2>>{log.err} 
mv $a.temp.vcf2 $a.temp.vcf
for chrid in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
    SnpSift -Xmx14g annotate -info AF,AF_nfe,grpmax,fafmax_faf95_max,nhomalt {params.gnomad}/gnomad.genomes.v4.0.sites.chr$chrid.vcf.bgz $a.temp.vcf > $a.temp.vcf2 2>>{log.err}
    mv $a.temp.vcf2 $a.temp.vcf
done
cat $a.temp.vcf|bgzip -c - > $b
tabix -p vcf $b
"""

rule Multiqc:
    input:
        expand(
            join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "{sample}.report.html"),
            sample=sampledic,
        )
    output:
        path   = directory(join(config['workdir'], "03.MultiQC")),
        report = join(config['workdir'], "03.MultiQC", "multiqc_report.html"),
    params:
        config = join(config['pipelinedir'], "multiqc.yaml"),
    log:
        out = join(config['pipelinedir'], "logs", "Multiqc", "Merge.o"),
        err = join(config['pipelinedir'], "logs", "Multiqc", "Merge.e"),
    threads:
        int(allocated("threads", "Multiqc", cluster))
    container:
        config['container']['multiqc']
    shell:
        '''
        cd {output.path}
        multiqc \
            -f \
            -c {params.config} \
            -o ./ \
            .. > {log.out} 2>{log.err}
        '''

rule Merge_summary:
    input:
        expand(
            join(config['workdir'], "02.Duplex", "{sample}", "{sample}.config.csv.summary.csv"),
            sample=sampledic,
        )
    output:
        path   = join(config['workdir'], "All.Duplex.summary.tsv"),
    params:
        path   = config['workdir']
    log:
        out = join(config['pipelinedir'], "logs", "Merge_summary", "Merge.o"),
        err = join(config['pipelinedir'], "logs", "Merge_summary", "Merge.e"),
    threads:
        int(allocated("threads", "Merge_summary", cluster))
    shell:
        '''
        cd {params.path}
        python {config[pipelinedir]}/scripts/merge_summary.py {input} > {output.path} 2>{log.err}
        '''
