#!/usr/bin/env python

import csv
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

   if sv_vcf.endswith('.gz'):
      print 'Input file must be an un-compressed, raw VCF file (extension \'.vcf\')'
      return
   f_in = open(sv_vcf,'r')
   f_out = open(out_tsv,'a')
   for line in f_in:
      if not line.startswith('#'):
         vcf_cols = line.rstrip().split('\t')
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

   f_in.close()
   f_out.close()


if __name__=="__main__": __main__()


