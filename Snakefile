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
        int(allocated("threads", "{rule}", cluster))
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
        int(allocated("threads", "{rule}", cluster))
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
        htmlout = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.QC.html"),
        jsonout = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.QC.json"),
    log:
        out = join(config['pipelinedir'], "logs", "Fastp", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Fastp", "{sample}.e"),
    threads:
        int(allocated("threads", "{rule}", cluster))
    container:
        config['container']['duplex']
    shell:
        '''
        fastp -i {input.reads1} \
            -I {input.reads2} \
            --stdout \
            -h {output.htmlout} \
            -j {output.jsonout} \
            -w {threads} \
            >> /dev/null 2> {log.err}
        '''

rule Duplex:
    input:
        reads1 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R1.fq.gz"),
        reads2 = join(config['workdir'], "01.MergeData", "{sample}", "{sample}.R2.fq.gz"),
    output:
        report = join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "{sample}.report.html"),
        summary= join(config['workdir'], "02.Duplex", "{sample}", "{sample}.config.csv.summary.csv"),
    params:
        path   = join(config['workdir'], "02.Duplex","{sample}","workdir"),
        config = join(config['pipelinedir'], "duplex_config", "{sample}.config.csv"),
    log:
        out = join(config['pipelinedir'], "logs", "Duplex", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Duplex", "{sample}.e"),
    threads:
        int(allocated("threads", "{rule}", cluster))
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

rule Multiqc:
    input:
        expand(
            join(config['workdir'], "02.Duplex", "{sample}", "workdir", "Final", "{sample}.report.html"),
            sample=sampledic,
        )
    output:
        path   = join(config['workdir'], "03.MultiQC"),
        report = join(config['workdir'], "03.MultiQC", "multiqc_report.html"),
    params:
        config = join(config['pipelinedir'], "multiqc.yaml"),
    log:
        out = join(config['pipelinedir'], "logs", "Multiqc", "Merge.o"),
        err = join(config['pipelinedir'], "logs", "Multiqc", "Merge.e"),
    threads:
        int(allocated("threads", "{rule}", cluster))
    container:
        config['container']['duplex']
    shell:
        cd {output.path}
        multiqc \
            -f \
            -c {params.config} \
            -o ./ \
            ..

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
        int(allocated("threads", "{rule}", cluster))
    shell:
        '''
        cd {params.path}
        cat {input}|head -1|sed "s/,/\t/g" > {output.path}
        cat {input}|grep -v "^RunID" |sed "s/,/\t/g">> {output.path}
        '''
