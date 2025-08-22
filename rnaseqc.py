#!/usr/bin/env python

# wrapper functions for filtering BAM files and QC plots for RNA-seq 
# Author: Yogesh Dhungana

import sys
import os
import argparse
from encode_lib_common import (
    copy_f_to_dir, log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext,
    strip_ext_bam, get_num_lines)
from encode_lib_genomic import (
    remove_chrs_from_bam, samtools_index,
    samtools_name_sort)
    

def parse_arguments():
    parser = argparse.ArgumentParser(prog='RNA-seq QC analysis and generation of QC plots.', description='')
    parser.add_argument('--bam', type=str, nargs=1, help='Path for raw BAM file.')
    parser.add_argument('--paired-end', action='store_true', help='Paired-end BAM.')
    parser.add_argument('--knowngeneBED', nargs=1, help='Path to Known gene bed for QC results.')
    parser.add_argument('--nth', type=int, default=8,help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='', type=str,help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    
    
    return args
    
    
    
def rnaseq_qc(bam, genebed, out_dir, nth):
    samtools_index(bam, nth, out_dir)
    basename = os.path.basename(strip_ext_bam(bam))
    prefix = os.path.join(out_dir, basename)
    bamstat = '{}.bamstat.txt'.format(prefix)
    readdist = '{}.readDist.txt'.format(prefix)
                          
    cmd1 = 'bam_stat.py -i {} > {}'
    cmd1 = cmd1.format(bam, bamstat)
    run_shell_cmd(cmd1)
    
    cmd2 = 'geneBody_coverage.py -r {} -i {} -o {}'
    cmd2 = cmd2.format(genebed, bam, prefix)
    run_shell_cmd(cmd2)
    
    cmd3 = 'inner_distance.py -r {} -i {} -o {}'
    cmd3 = cmd3.format(genebed, bam, prefix)
    run_shell_cmd(cmd3)
    
    cmd4 = 'junction_annotation.py -r {} -i {} -o {}'
    cmd4 = cmd4.format(genebed, bam, prefix)
    run_shell_cmd(cmd4)
    
    cmd5 = 'junction_saturation.py -r {} -i {} -o {}'
    cmd5 = cmd5.format(genebed, bam, prefix)
    run_shell_cmd(cmd5)
    
    cmd6 = 'read_distribution.py -r {} -i {} > {}'
    cmd6 = cmd6.format(genebed, bam, readdist)
    run_shell_cmd(cmd6)
    
    cmd7 = 'read_duplication.py -i {} -o {}'
    cmd7 = cmd7.format(bam, prefix)
    run_shell_cmd(cmd7)
    
    cmd8 = 'read_GC.py -i {} -o {}'
    cmd8 = cmd8.format(bam, prefix)
    run_shell_cmd(cmd8)
    
    cmd9 = 'read_NVC.py -i {} -o {}'
    cmd9 = cmd9.format(bam, prefix)
    run_shell_cmd(cmd9)
    
    cmd10 = 'read_quality.py -i {} -o {}'
    cmd10 = cmd10.format(bam, prefix)
    run_shell_cmd(cmd10)
    
    
def main():
    args = parse_arguments()
    log.info('Initializing and creating Output directory...')
    mkdir_p(args.out_dir)
    log.info("Generating QC plots")
    rnaseq_qc(args.bam[0], args.knowngeneBED[0], args.out_dir, args.nth)
    log.info("All Done.")
        

if __name__ == '__main__':
    main()