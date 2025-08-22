#!/usr/bin/env python

# Script to submit RNA-seq jobs in parallel
# Author Yogesh Dhungana

import glob
import subprocess
import re
import os
import sys
import argparse
sys.path.append("/home/ydhungan/pipelines/rnaseq_pipeline/bin")
from encode_lib_common import (
    copy_f_to_dir, copy_f_to_f, log, ls_l, mkdir_p, read_tsv, rm_f,
    run_shell_cmd, strip_ext_fastq)

def parse_arguments():
    parser = argparse.ArgumentParser(prog='Wrapper for BSUB job submission for RNA-seq data.', description='')
    parser.add_argument('--path-to-fastqs', default='', type=str, help= 'Path to FASTQ files.')
    parser.add_argument('--paired-end', action='store_true', help='Paired-end fastqs.')
    parser.add_argument('--configfile',default='/home/ydhungan/pipelines/rnaseq_pipeline/bin/config.ini', type=str, help='Configuration file with information for genome database')
    parser.add_argument('--adapters',default='/home/ydhungan/adapters/TruSeqAdapters.fa', type=str, help='Fasta file of adapter sequences')
    parser.add_argument('--genome', default='mm10', choices=['mm10','hg19','hg38'], help='Genome to use for quantification')
    parser.add_argument('--method',default='salmon', choices=['salmon','rsem'], help='method of RNA-seq quantification')
    parser.add_argument('--strandedness',default='none', choices=['none','forward','reverse'], help='method of RNA-seq quantification')
    parser.add_argument('--trim-reads',default=1,type=int,choices=[1,0], help='Option to trim or not to trim the input reads')
    parser.add_argument('--memory', default='16000', type=str, help= 'Memory requested to run the analysis.')
    parser.add_argument('--queue', default='standard', type=str,help='Queue to submit the job in HPCF (use bqueues to choose).')
    parser.add_argument('--out-dir', default='.', type=str,help='Output Directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    
    return args
    
args = parse_arguments()    
log.info("Making output directory...")
mkdir_p(args.out_dir)

if args.paired_end: 
    fastqR1 = glob.glob(args.path_to_fastqs+'/*R1_001.fastq.gz')
    fastqR2 = glob.glob(args.path_to_fastqs+'/*R2_001.fastq.gz')
    fastqR1.sort()
    fastqR2.sort()
else: 
    fastqR1 = glob.glob(args.path_to_fastqs+'/*R1_001.fastq.gz')
    fastqR1.sort()


hpcfsubmit = 'bsub ' + '-R ' + '"rusage[mem=' + args.memory + ']" ' + '-q ' + args.queue + ' < '
#print 'Job submitted with command: {}'.format(hpcfsubmit)
def submit_job(jobf):
  os.system('{}'.format(hpcfsubmit) + jobf)
  

def create_job_file_pe(samplefile1, samplefile2, adapters, genome, method, strandedness, out_dir, trim_reads):
    basename1 = os.path.basename(strip_ext_fastq(samplefile1))
    basename2 = os.path.basename(strip_ext_fastq(samplefile2))
    prefix = os.path.join(out_dir,basename1)
    job_header = '#!/bin/bash\n'
    job_header += '#BSUB -P RNAseq\n'
    job_header += '#BSUB -J {}_RNAseq\n'
    job_header += '#BSUB -oo {}'+'/rnaseqlog.out\n'
    job_header += '#BSUB -eo {}'+'/rnaseqerror.err\n'
    job_header += '#BSUB -n 1\n'
    job_header += '#BSUB -N ydhungan@stjude.org\n'
    job_header = job_header.format(basename1,prefix, prefix)
    
    
    ### Load all the required module for analysis: 
    module1 = 'module load trimmomatic/0.36\n'
    module1 += 'module load samtools/1.9\n'
    module1 += 'module load rsem/1.3.1\n'
    module1 += 'module load python/3.7.0\n'
    job_body1 = ''
    job_body2 = ''
    ### Write job body to run each wrapper for the sample:
    if trim_reads==1:
        job_body1 = 'python /home/ydhungan/pipelines/rnaseq_pipeline/bin/trimmomatic_run.py --fastqs {} {} --nth 8 --adapters {} --out-dir {} --paired-end --trim-reads {}'
        job_body1 = job_body1.format(samplefile1, samplefile2, adapters, prefix, trim_reads)
        job_body2 = 'python /home/ydhungan/pipelines/rnaseq_pipeline/bin/rnaquantification.py '
        job_body2 += '--fastqs {} {} --genome {} --method {} --strandedness {} --paired-end --out-dir {}'
        job_body2 = job_body2.format(prefix+'/'+basename1 +'.trim.paired.fastq.gz', prefix+'/'+basename2 +'.trim.paired.fastq.gz', genome, method, strandedness, prefix)
    else:
        job_body1 = 'python /home/ydhungan/pipelines/rnaseq_pipeline/bin/trimmomatic_run.py --fastqs {} {} --nth 8 --adapters {} --out-dir {} --paired-end --trim-reads {}'
        job_body1 = job_body1.format(samplefile1, samplefile2, adapters, prefix, trim_reads)
        job_body2 = 'python /home/ydhungan/pipelines/rnaseq_pipeline/bin/rnaquantification.py '
        job_body2 += '--fastqs {} {} --genome {} --method {} --strandedness {} --paired-end --out-dir {}'
        job_body2 = job_body2.format(samplefile1, samplefile2, genome, method, strandedness, prefix)
    
    jobfile = prefix+".sh"
    with open(jobfile,"w") as new_file:
        new_file.write(job_header+module1+'\n'+job_body1+'\n'+job_body2+'\n')
    
    return jobfile
    
def create_job_file_se(samplefile1, adapters, genome, method, strandedness, out_dir, trim_reads):
    basename = os.path.basename(strip_ext_fastq(samplefile1))
    samplename = basename.split('_')[0]
    prefix = os.path.join(out_dir,basename)
    job_header = '#!/bin/bash\n'
    job_header += '#BSUB -P RNAseq\n'
    job_header += '#BSUB -J {}_RNAseq\n'
    job_header += '#BSUB -oo {}'+'/rnaseqlog.out\n'
    job_header += '#BSUB -eo {}'+'/rnaseqerror.err\n'
    job_header += '#BSUB -n 1\n'
    job_header += '#BSUB -N ydhungan@stjude.org\n'
    job_header = job_header.format(basename,prefix, prefix)
    
    
    ### Load all the required module for analysis: 
    module1 = 'module load trimmomatic/0.36\n'
    module1 += 'module load samtools/1.9\n'
    module1 += 'module load rsem/1.3.1\n'
    module1 += 'module load python/3.7.0\n'
    job_body1 = ''
    job_body2 = ''
    ### Write job body to run each wrapper for the sample:
    if trim_reads==1:
        job_body1 = 'python /home/ydhungan/pipelines/rnaseq_pipeline/bin/trimmomatic_run.py --fastqs {} --nth 8 --adapters {} --out-dir {} --trim-reads {}'
        job_body1 = job_body1.format(samplefile1, adapters, prefix, trim_reads)
        job_body2 = 'python /home/ydhungan/pipelines/rnaseq_pipeline/bin/rnaquantification.py '
        job_body2 += '--fastqs {} --genome {} --method {} --strandedness {} --out-dir {}'
        job_body2 = job_body2.format(prefix+'/'+samplename +'_R1_001.trim.paired.fastq.gz', genome, method, strandedness, prefix)
    else:
        job_body1 = 'python /home/ydhungan/pipelines/rnaseq_pipeline/bin/trimmomatic_run.py --fastqs {} --nth 8 --adapters {} --out-dir {} --trim-reads {}'
        job_body1 = job_body1.format(samplefile1, adapters, prefix, trim_reads)
        job_body2 = 'python /home/ydhungan/pipelines/rnaseq_pipeline/bin/rnaquantification.py '
        job_body2 += '--fastqs {} --genome {} --method {} --strandedness {} --out-dir {}'
        job_body2 = job_body2.format(samplefile1, genome, method, strandedness, prefix)
    
    jobfile = prefix+".sh"
    with open(jobfile,"w") as new_file:
        new_file.write(job_header+module1+'\n'+job_body1+'\n'+job_body2+'\n')
    
    return jobfile



if args.paired_end: 
    for fastq in range(0,len(fastqR1)):
        submit_job(create_job_file_pe(fastqR1[fastq], fastqR2[fastq],args.adapters, args.genome, args.method, args.strandedness, args.out_dir, args.trim_reads))
        
else: 
    for fastq in range(0,len(fastqR1)):
        submit_job(create_job_file_se(fastqR1[fastq], args.adapters, args.genome, args.method, args.strandedness, args.out_dir, args.trim_reads))
    
    
    