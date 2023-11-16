#!/usr/bin/env python
import os, sys
import gzip
import yaml
import time
import pandas as pd
import subprocess
from optparse import OptionParser

def main():
    parser = OptionParser(usage="usage: %prog [options]")
    parser.add_option("-s","--samplesheet", dest="samplesheet",default=None, help="Samplesheet to be used for the pipeline.")
    parser.add_option("-w","--workdir",dest="workdir",default=None, help="Directory to run the pipeline.")
    parser.add_option("-c","--cachedir",dest="cachedir",default=None,help="Directory of the cache to be used. Default: /Path/to/pipeline/cache")
    parser.add_option("-C","--cacheinit", dest="cacheinit", action = "store_true", default=False, help="Run the cache initiation.")
    parser.add_option("-u","--unlock", dest="unlock", action = "store_true", default=False, help="Unlock the workdir from previous abnormal pipeline.")
    parser.add_option("-S","--silent", dest="silent", action = "store_true", default=False, help="Keep silent, stop sending Slurm email notice.")
    parser.add_option("-D","--dry", dest="dryrun", action = "store_true", default=False, help="Perform snakemake dryrun.")
    parser.add_option("-R","--region",dest="targetbed",default=None,help="Directory to the targeted region bed file.") 
    parser.add_option("-E","--exclude",dest="maskbed",default=None,help="Directory to the exclude region bed file.")
    parser.add_option("-t","--twinstrand", action = "store_true",dest="chemistry_ts",default=None,help="Use Twinstrand chemistry. [Default]")
    parser.add_option("-T","--twist", action = "store_true",dest="chemistry_tw",default=None,help="Use Twist chemistry.")
    parser.add_option("-U","--umilen",dest="umilen",default=8,help="UMI length. [Default=8]")
    parser.add_option("-K","--spacerlen",dest="spacerlen",default=2,help="Spacer length. [Default=2]")
    parser.add_option("-L","--loclen",dest="loclen",default=0,help="Location identifier length. [Default=0]")
    parser.add_option("-B","--clipbegin",dest="clipbegin",default=0,help="Hard clip first N bases. [Default=0]")
    parser.add_option("-b","--clipend",dest="clipend",default=0,help="Hard clip last N bases. [Default=0]")
    
    (options, args) = parser.parse_args()

    if not options.chemistry_ts and not options.chemistry_tw:
        options.chemistry_ts = True
    if options.chemistry_ts:
        chemistry = 'ts'
        adapterseq = 'NNNNNNNNNNAGATCGGAAGAG'
        if options.chemistry_tw:
            raise Exception("--twinstrand and --twist cannot be specified together.")
    elif options.chemistry_tw:
        chemistry = 'tw'
        adapterseq = 'NNNNNNNAGATCGGAAGAG'
        options.umilen=5
        options.spacerlen=2
        
    

    #######################################################################################
    # Cache mode
    #######################################################################################
    ## Set default cachedir inside the pipeline folder
    pipedir = os.path.dirname(__file__)
    if not options.cachedir:
        options.cachedir = os.path.join(pipedir,'cache')
    elif options.cachedir[0] != '/':
        print('Provided cachedir is not an absolute path. Placing it inside the pipeline dir:')
        print(os.path.join(pipedir, options.cachedir))
        options.cachedir = os.path.join(pipedir, options.cachedir)
        
    ## Do cache if cacheinit provided
    if options.cacheinit:
        print('Cache initiation mode, ignore all options except -c (--cachedir)')
        print('Init cache to %s...'%(options.cachedir))
        from scripts.cache import container,reference
        
        # Init container
        with open(os.path.join(pipedir, 'config', 'container.yaml'), 'r') as infile:
            config = yaml.safe_load(infile)
        config['pipelinedir'] = pipedir
        config['cachedir'] = options.cachedir
        status_cache_1 = container(config)
        if not status_cache_1:
            raise Exception("Container cache failed.")
        
        # Init reference
        with open(os.path.join(pipedir, 'config', 'reference.yaml'), 'r') as infile:
            config = yaml.safe_load(infile)
        config['pipelinedir'] = pipedir
        config['cachedir'] = options.cachedir
        status_cache_2 = reference(config)
        if not status_cache_2:
            raise Exception("Reference cache failed.")
        return
    
    #######################################################################################
    # Unlock mode
    #######################################################################################
    if options.unlock:
        print('Pipeline unlock mode, ignore all options except -w (--workdir)')
        
        if os.path.exists(os.path.join(options.workdir,'.snakemake','locks','0.input.lock')):
            os.remove(os.path.join(options.workdir,'.snakemake','locks','0.input.lock'))
        if os.path.exists(os.path.join(options.workdir,'.snakemake','locks','0.output.lock')):
            os.remove(os.path.join(options.workdir,'.snakemake','locks','0.output.lock'))
            
        print('Pipeline directory:',pipedir, 'unlocked!')
        return

    #######################################################################################
    # Prepare workdir
    #######################################################################################
    if not options.targetbed:
        raise Exception("--region must be specified.")
    
    snapshot = 'run_%s'%(time.strftime("%Y%m%d%H%M%S"))
    from scripts.load import samplesheet
    if not options.samplesheet:
        parser.error("Samplesheet (-s) is not specified.")
    if not options.workdir:
        parser.error("Workdir (-w) is not specified.")
    if not os.path.exists(os.path.join(options.cachedir, 'reference', 'ref.ok')):
        parser.error("Cache (%s) is not ready for use, check parameters!"%(options.cachedir))
    if not os.path.exists(os.path.join(os.path.realpath(options.workdir), 'Pipe_runtime', snapshot, 'logs', 'slurm')):
        os.system('mkdir -p %s'%(os.path.join(os.path.realpath(options.workdir), 'Pipe_runtime', snapshot, 'logs', 'slurm')))
    # Check and write samplesheet to the workdir
    try:
        sample, run, run2sample = samplesheet(options.samplesheet)
    except:
        parser.error("Samplesheet failed to load. Check %s"%(options.samplesheet))
    
    #######################################################################################
    # Prepare config file
    #######################################################################################
    config = {}
    config['samplesheet'] = os.path.realpath(options.samplesheet)
    config['workdir'] = os.path.realpath(options.workdir)
    config['snapshot'] = snapshot
    config['pipelinedir'] = os.path.join(config['workdir'], 'Pipe_runtime', snapshot)
    config['cachedir'] = os.path.realpath(options.cachedir)
    config['options'] = {}
               
    for cfg in ['config.yaml','reference.yaml']:
        with open(os.path.join(pipedir, 'config', cfg), 'r') as infile:
            tmpcfg = yaml.safe_load(infile)
            for i in tmpcfg:
                if i not in config:
                    config[i] = tmpcfg[i]
                    
    with open(os.path.join(pipedir, 'config', 'container.yaml'), 'r') as infile:
        tmpcfg = yaml.safe_load(infile)
        config['container'] = {}
        for c in tmpcfg['docker']:
            config['container'][c] = os.path.join(config['cachedir'], 'container', '%s.simg'%(c))
        for c in tmpcfg['simg']:
            config['container'][c] = os.path.join(config['cachedir'], 'container', '%s.simg'%(c))
     
    with open(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'config.yaml'), 'w') as savefile:
        savefile.write(yaml.dump(config))
    os.system('cp %s %s'%(options.samplesheet, os.path.join(options.workdir, 'Pipe_runtime', snapshot, 'samplesheet.tsv')))
    
    # copy cluster.yaml
    if os.path.exists(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'cluster.yaml')):
        os.remove(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'cluster.yaml'))
    os.system('cp %s %s'%(os.path.join(pipedir, 'config', 'cluster.yaml'), os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'cluster.yaml')))
    
    #######################################################################################
    # Prepare single duplex-seq config file
    ####################################################################################### 
    from scripts.utils import get_read_length
    for s in sample:
        ssdf = pd.read_csv(os.path.join(pipedir, 'duplex.config.csv'))
        os.system('mkdir -p %s'%(os.path.join(options.workdir, 'Pipe_runtime', snapshot, 'duplex_config')))
        ssdf.loc[s,'sample'] = s
        ssdf.loc[s,'rglb'] = s
        ssdf.loc[s,'rgpl'] = 'illumina'
        ssdf.loc[s,'rgpu'] = s
        ssdf.loc[s,'rgsm'] = s
        ssdf.loc[s,'reference'] = os.path.join(config['cachedir'],'reference','hg38','Homo_sapiens_assembly38.fasta')
        ssdf.loc[s,'target_bed'] = os.path.realpath(options.targetbed)
        ssdf.loc[s,'maskBed'] = 'NONE'
        ssdf.loc[s,'blast_db'] = 'none'
        ssdf.loc[s,'targetTaxonId'] = '9606'
        ssdf.loc[s,'baseDir'] = 'workdir'
        ssdf.loc[s,'in1'] = 'read1.fq.gz'
        ssdf.loc[s,'in2'] = 'read2.fq.gz'
        ssdf.loc[s,'mqFilt'] = 0
        ssdf.loc[s,'minMem'] = 3
        ssdf.loc[s,'maxMem'] = 200
        ssdf.loc[s,'cutOff'] = 0.7
        ssdf.loc[s,'nCutOff'] = 0.02
        ssdf.loc[s,'umiLen'] = options.umilen
        ssdf.loc[s,'spacerLen'] = options.spacerlen
        ssdf.loc[s,'locLen'] = options.loclen
        ssdf.loc[s,'readLen'] = get_read_length(run[sample[s][0]]['r1'])
        ssdf.loc[s,'adapterSeq'] = adapterseq 
        ssdf.loc[s,'clipBegin'] = options.clipbegin
        ssdf.loc[s,'clipEnd'] = options.clipend
        ssdf.loc[s,'minClonal'] = 0
        ssdf.loc[s,'maxClonal'] = 0.1
        ssdf.loc[s,'minDepth'] = 100
        ssdf.loc[s,'maxNs'] = 1
        ssdf.loc[s,'recovery'] = 'noRecovery.sh'
        ssdf.loc[s,'cluster_dist'] = 10
        ssdf.loc[s,'cm_outputs'] = 'GB'
        ssdf.loc[s,'cm_sumTypes'] = 'GT'
        ssdf.loc[s,'cm_filters'] = 'none'
        ssdf.loc[s,'runSSCS'] = 'TRUE'
        ssdf.loc[s,'rerun_type'] = 0
        ssdf.to_csv(os.path.join(options.workdir, 'Pipe_runtime', snapshot, 'duplex_config','%s.config.csv'%s), index=None)

    #######################################################################################
    # Prepare other files
    #######################################################################################
    
    # copy multiqc config
    #if os.path.exists(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'multiqc.yaml')):
    #    os.remove(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'multiqc.yaml'))
    #os.system('cp %s %s'%(os.path.join(pipedir, 'config', 'multiqc.yaml'), os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'multiqc.yaml')))
    
    # copy scripts
    if os.path.exists(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'scripts')):
        os.system('rm -rf %s'%(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'scripts')))
    os.system('cp -r %s %s'%( os.path.join(pipedir, 'scripts'), os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'scripts') ))
    
    # copy snakefile
    if os.path.exists(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'Snakefile')):
        os.system('rm %s'%(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'Snakefile')))
    os.system('cp -r %s %s'%( os.path.join(pipedir, 'Snakefile'), os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'Snakefile') ))
    
    # copy multiqc config
    if os.path.exists(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'multiqc.yaml')):
        os.system('rm %s'%(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'multiqc.yaml')))
    os.system('cp -r %s %s'%( os.path.join(pipedir, "config", 'multiqc.yaml'), os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'multiqc.yaml') ))
    
    # write pipeline submission bash script: captain.sh
    with open(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'captain.sh'), 'w') as savefile:
        with open(os.path.join(pipedir, 'captain.sh')) as infile:
            captain = infile.read()
            if options.silent:
                captain = captain.replace('#SBATCH --mail-type=BEGIN,END,FAIL\n','')
            captain = captain.replace('[[PIPENICKNAME]]','Duplex_seq')
            captain = captain.replace('[[WORKDIR]]', config['workdir'])
            captain = captain.replace('[[SNAPSHOT]]', config['snapshot'])
            captain = captain.replace('[[BINDPATH]]', config['bindpath'])
            savefile.write(captain)
    os.chmod(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'captain.sh'), 0o755)
    
    #######################################################################################
    # submit captain.sh to slurm
    #######################################################################################
    if os.path.exists(os.path.join(config['workdir'], '.snakemake', 'locks', '0.input.lock')) or os.path.exists(os.path.join(config['workdir'], '.snakemake', 'locks', '0.output.lock')):
        raise Exception("Workdir locked! Run the pipemaster with -u/--unlock to unlock the workdir first.")

    if options.dryrun:
        os.system('cd %s && snakemake -n'%(config['pipelinedir']))
    else:
        cmd = 'cd %s && '%config['pipelinedir']
        cmd += 'sbatch captain.sh'
        try:
            subprocess.check_call([cmd], shell=True, executable="/bin/bash")
            print('Pipeline submitted!')
        except:
            raise Exception("Pipeline submission failed.") 

if __name__ == '__main__':
    main()