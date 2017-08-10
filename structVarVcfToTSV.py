#!/usr/bin/env python

import re
import argparse
import os
import logging
import sys

def __main__():
   
   parser = argparse.ArgumentParser(description='Transform VCF output from structural variant calling (i.e. Manta) to TSV format accepted by Variant Effect Predictor (VEP)')
   parser.add_argument('sv_vcf', help='Structural variant VCF input file with somatic query variants')
   parser.add_argument('out_tsv', help='Structural variants in VEP-compatible input format (tab-separated values)')
   args = parser.parse_args()
   
   sv_vcf_VEP_compatible_input(args.sv_vcf, args.out_tsv)

def sv_vcf_VEP_compatible_input(sv_vcf, out_tsv):

   logger = getlogger('sv_vcf_transformation')
   if sv_vcf.endswith('.gz'):
      logger.info('Input file must be an un-compressed, raw VCF file (extension \'.vcf\')')
      return
   
   logger.info('Reading input file: ' + str(sv_vcf))
   f_in = open(sv_vcf,'r')
   f_out = open(out_tsv,'a')
   n_pass = 0
   n_non_pass = 0
   for line in f_in:
      if not line.startswith('#'):
         vcf_cols = line.rstrip().split('\t')
         if len(vcf_cols) >= 8:
            vcf_info_data = vcf_cols[7]
            vcf_filter = vcf_cols[6]
            vcf_alt = vcf_cols[4]
            if re.search(r'SVTYPE=(INS|DEL|DUP|TDUP)',vcf_info_data):
               sv_entry_chrom = str(vcf_cols[0])
               sv_entry_start = str(vcf_cols[1])
               sv_entry_stop = -1
               sv_entry_type = '_na'
               
               vcf_info_tags = vcf_info_data.split(';')
               for t in vcf_info_tags:
                  if t.startswith('SVTYPE='):
                     sv_entry_type = str(t.split('=')[1])
                     
                  
                  if t.startswith('END='):
                     sv_entry_stop = str(t.split('=')[1])
            
               if sv_entry_type == 'DUP' and vcf_alt == '<DUP:TANDEM>':
                  sv_entry_type = 'TDUP'
               
               if vcf_filter == 'PASS':    
                  f_out.write(str(sv_entry_chrom) + '\t' + str(sv_entry_start) + '\t' + str(sv_entry_stop) + '\t' + str(sv_entry_type) + '\n')
                  n_pass += 1
               else:
                  n_non_pass += 1
               
   logger.info('Number of structural variant calls - PASS: ' + str(n_pass))
   logger.info('Number of structural variant calls - non-PASS: ' + str(n_non_pass))
   logger.info('VEP-compatible, tab-separated file in: ' + str(out_tsv))

   f_in.close()
   f_out.close()

def getlogger(logger_name):
	logger = logging.getLogger(logger_name)
	logger.setLevel(logging.DEBUG)

	# create console handler and set level to debug
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.DEBUG)

	# add ch to logger
	logger.addHandler(ch)
	
	# create formatter
	formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", "20%y-%m-%d %H:%M:%S")
	
	#add formatter to ch
	ch.setFormatter(formatter)
	
	return logger

if __name__=="__main__": __main__()


