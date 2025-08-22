#!/usr/bin/env python

# adapter trimming wrapper for all seq experiments using Trimmomatic v 0.36
# Author: Yogesh Dhungana

import sys
import os 
import argparse
import copy 
from encode_lib_common import (
    copy_f_to_dir, copy_f_to_f, log, ls_l, mkdir_p, read_tsv, rm_f,
    run_shell_cmd, strip_ext_fastq)

def parse_arguments(debug=False):
    parser = argparse.ArgumentParser(prog = 'Adpater trimmer', description='')
    parser.add_argument('--fastqs', type=str,nargs='*', help='Reads of FASTQs must be compressed with gzip (.gz)')
    parser.add_argument('--trimmomatic-parm', type=str, default='LEADING:10 TRAILING:10 SLIDINGWINDOW:4:18 MINLEN:25', help= 'default: LEADING:10 TRAILING:10 SLIDINGWINDOW:4:18 MINLEN:25')
    parser.add_argument('--nth', type=str, default=1, help='Number of threads to parallelize.')
    parser.add_argument('--adapters',nargs=1, type=str, help='Fasta file with adapter sequences (.fa)')
    parser.add_argument('--out-dir', default='', type=str, help='Output Directory')
    parser.add_argument('--paired-end', action='store_true', help='Paired-end FASTQs')
    parser.add_argument('--trim-reads', default=1,type=int, choices=[1,0], help='Choose to trim the adapter sequence')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    # Parse fastq commandline for paired end and single end: 
    if args.paired_end and len(args.fastqs)!=2:
        raise argparse.ArgumentTypeError(
            'Need 2 fastq files for paired end read')
    if not args.paired_end and len(args.fastqs)!=1:
        raise argparse.ArgumentTypeError(
            'Need 1 fastq file for single end read')
    return (args)


def trim_adapter_se(fastq, adapter, nthread, trimmomatic_parm, out_dir, trim_reads):
    prefix = os.path.join(out_dir,
                    os.path.basename(strip_ext_fastq(fastq)))
    if trim_reads==1:
        trimmed1 = '{}.trim.paired.fastq.gz'.format(prefix)
        cmd = 'trimmomatic SE -phred33 -threads ' + nthread + ' ' + fastq + ' ' + trimmed1 + ' ILLUMINACLIP:' + adapter + ':2:30:10 ' + trimmomatic_parm 
        run_shell_cmd(cmd)
    else:
        mkdir_p(prefix)
        copy_f_to_dir(fastq,prefix)

def trim_adapter_pe(fastq1, fastq2, adapter, nthread,trimmomatic_parm, out_dir, trim_reads):
    prefix = os.path.join(out_dir,
                    os.path.basename(strip_ext_fastq(fastq1)))
    if trim_reads==1:
        prefix1 = os.path.join(out_dir, os.path.basename(strip_ext_fastq(fastq1)))
        prefix2 = os.path.join(out_dir, os.path.basename(strip_ext_fastq(fastq2)))
        trimmed1 = '{}.trim.paired.fastq.gz'.format(prefix1)
        trimmed2 = '{}.trim.unpaired.fastq.gz'.format(prefix1)
        trimmed3 = '{}.trim.paired.fastq.gz'.format(prefix2)
        trimmed4 = '{}.trim.unpaired.fastq.gz'.format(prefix2)
        cmd = 'trimmomatic PE -phred33 -threads ' + nthread + ' ' + fastq1 + ' '+ fastq2 + ' ' + trimmed1 + ' ' + trimmed2 +' ' + trimmed3 + ' ' + trimmed4 +  ' ILLUMINACLIP:' + adapter + ':2:30:10 ' + trimmomatic_parm
        run_shell_cmd(cmd)
    else:
        mkdir_p(prefix)
        copy_f_to_dir(fastq1, prefix)
        copy_f_to_dir(fastq2, prefix)


def main():
    # read param
    args = parse_arguments()
    log.info('Initializing and creating Output directory...')
    mkdir_p(args.out_dir)
    fastqs = args.fastqs
    adapters = args.adapters
    trimmomatic_parm = args.trimmomatic_parm
    nthread = args.nth
    out_directory = args.out_dir
    if args.paired_end:
        log.info('Trimming reads in Paired end mode...')
        trim_adapter_pe(fastqs[0], fastqs[1], adapters[0], nthread, trimmomatic_parm, out_directory, args.trim_reads)
    else:
        log.info('Trimming reads in Single end mode...')
        trim_adapter_se(fastqs[0], adapters[0], nthread, trimmomatic_parm, out_directory, args.trim_reads)
    log.info('All done.')

if __name__ == '__main__':
    main() 
    
    

