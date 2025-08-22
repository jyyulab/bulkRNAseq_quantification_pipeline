#!/usr/bin/env python

# Quantification wrapper for RNA-seq experiments
# Author: Yogesh Dhungana

import sys
import os 
import argparse
import copy 
import configparser
from encode_lib_common import (
    copy_f_to_dir, copy_f_to_f, log, ls_l, mkdir_p, read_tsv, rm_f,
    run_shell_cmd, strip_ext_fastq)
    

def parse_arguments(debug=False):
    parser = argparse.ArgumentParser(prog = 'RNA-seq Quantification using salmon or RSEM', description='')
    parser.add_argument('--fastqs', type=str, nargs='*', help='Reads of FASTQs must be compressed with gzip (.gz)')
    parser.add_argument('--genome', default='mm10', choices=['mm10','hg19','hg38'], help='Genome to use for quantification')
    parser.add_argument('--method',default='salmon', choices=['salmon','rsem'], help='method of RNA-seq quantification')
    parser.add_argument('--strandedness',default='none', choices=['none','forward','reverse'], help='method of RNA-seq quantification')
    parser.add_argument('--configfile',default='/home/ydhungan/pipelines/rnaseq_pipeline/bin/config.ini', type=str, help='Configuration file with information for genome database')
    parser.add_argument('--nth', type=str, default=8, help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='', type=str, help='Output Directory')
    parser.add_argument('--paired-end', action='store_true', help='Paired-end FASTQs')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    # Parse fastq commandline for paired end and single end: 
    if args.paired_end and len(args.fastqs) < 2:
        raise argparse.ArgumentTypeError(
            'Need 2 fastq files for paired end read')
    if not args.paired_end and len(args.fastqs) < 1:
        raise argparse.ArgumentTypeError(
            'Need 1 fastq file for single end read')
    return args
    
def parseConfig(configfile, genome, method):
    Config = configparser.ConfigParser()
    Config.read(configfile)
    # when the method is RSEM
    if genome == 'mm10' and method == 'rsem':
        path_to_genome = Config.get('mm10' , 'rsem')
        path_to_annotation = ''
        tool_path = Config.get('mm10' , 'star_path')
    elif genome == 'hg19' and method == 'rsem':
        path_to_genome = Config.get('hg19' , 'rsem')
        path_to_annotation = ''
        tool_path = Config.get('hg19' , 'star_path')
    elif genome == 'hg38' and method == 'rsem':
        path_to_genome = Config.get('hg38' , 'rsem')
        path_to_annotation = ''
        tool_path = Config.get('hg38' , 'star_path')
    # when method is salmon
    elif genome == 'mm10' and method == 'salmon':
        path_to_genome = Config.get('mm10' , 'salmon_index')
        path_to_annotation = Config.get('mm10' , 'salmon_annotation')
        tool_path = Config.get('mm10' , 'salmon_path')
    elif genome == 'hg19' and method == 'salmon':
        path_to_genome = Config.get('hg19' , 'salmon_index')
        path_to_annotation = Config.get('hg19' , 'salmon_annotation')
        tool_path = Config.get('hg19' , 'salmon_path')
    elif genome == 'hg38' and method == 'salmon':
        path_to_genome = Config.get('hg38' , 'salmon_index')
        path_to_annotation = Config.get('hg38' , 'salmon_annotation')
        tool_path = Config.get('hg38' , 'salmon_path')

    return (path_to_genome, path_to_annotation, tool_path)
    
    
def rsem_pe(fastq1, fastq2, nth, path_to_genome, tool_path, strandedness, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq1))
    prefix = os.path.join(out_dir, basename)
    cmd = 'rsem-calculate-expression --num-threads {} --star --star-path {} --strandedness {} --sort-bam-by-coordinate --phred33-quals --paired-end {} {} {} {}'
    cmd = cmd.format(nth, tool_path, strandedness, fastq1, fastq2, path_to_genome, prefix)
    run_shell_cmd(cmd)
    
def rsem_se(fastq, nth, path_to_genome, tool_path, strandedness, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir, basename)
    cmd = 'rsem-calculate-expression --num-threads {} --star --star-path {} --strandedness {} --sort-bam-by-coordinate --phred33-quals {} {} {}'
    cmd = cmd.format(nth, tool_path, strandedness, fastq, path_to_genome, prefix)
    run_shell_cmd(cmd)    
    

def salmon_pe(fastq1, fastq2, nth, path_to_genome,path_to_annotation, tool_path, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq1))
    prefix = os.path.join(out_dir, basename)
    cmd = tool_path
    cmd += ' quant -i {} -l A -p {} -g {} -1 {} -2 {} -o {}'
    cmd = cmd.format(path_to_genome, nth, path_to_annotation, fastq1, fastq2, prefix)
    run_shell_cmd(cmd)
def salmon_se(fastq, nth, path_to_genome, path_to_annotation, tool_path, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir, basename)
    mkdir_p(prefix)
    cmd = tool_path
    cmd += ' quant -i {} -l A -p {} -g {} -r {} -o {}'
    cmd = cmd.format(path_to_genome, nth, path_to_annotation, fastq, prefix)
    run_shell_cmd(cmd)
    
def main():
    args = parse_arguments()
    log.info('Initializing and creating Output directory...')
    mkdir_p(args.out_dir)
    log.info('Parsing Config file for genome location and prebuild indexes...')
    path_to_genome, path_to_annotation, tool_path = parseConfig(args.configfile, args.genome, args.method)
    if args.paired_end and args.method=='salmon':
        log.info('Quantification using salmon in paired end mode...')
        salmon_pe(args.fastqs[0], args.fastqs[1], args.nth, path_to_genome, path_to_annotation, tool_path, args.out_dir)
    if not args.paired_end and args.method=='salmon':
        log.info('Quantification using salmon in single end mode...')
        salmon_se(args.fastqs[0], args.nth, path_to_genome, path_to_annotation, tool_path, args.out_dir)
    if args.paired_end and args.method=='rsem':
        log.info('Quantification using RSEM in paired end mode...')
        rsem_pe(args.fastqs[0], args.fastqs[1], args.nth, path_to_genome, tool_path, args.strandedness, args.out_dir)
    if not args.paired_end and args.method=='rsem':
        log.info('Quantification using RSEM in single end mode...')
        rsem_se(args.fastqs[0], args.nth, path_to_genome, tool_path, args.strandedness, args.out_dir)
    log.info("All Done.")

if __name__ == '__main__':
    main()
    