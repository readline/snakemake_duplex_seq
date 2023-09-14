import os
import pandas as pd
from os.path import join
from scripts.load import samplesheet

snakedir = os.getcwd()
print(snakedir)
configfile: 'config.yaml'
for s in config['simg']:
    config['simg'][s] = snakedir+'/'+config['simg'][s]
print(config)
sampledic, libdic, rundic = samplesheet(config['samplesheet'])
workdir: config['workdir']

rule all:
    input:

include: join("rules", "qc.smk")