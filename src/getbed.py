'''
A tool extract gene/burden from ExAC vcf files
Author: Khalid Mahmood
Contact: khalid.mahmood@unimelb.edu.au
Copyright: 2015
'''

#!/usr/bin/python

from annotations import getTabixValCondel
import sys
import os
import argparse
import getopt
import vcf
import re
import array
import pysam

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--vcf", type=str, dest="vcf", help="Input variant file (vcf)", required=True)
    parser.add_argument("-o", "--output", type=str, dest="out", help="Output file (tabular)", required=True)
    #parser.add_argument("-X", "--exac", type=str, dest="exac_af_threshold", help="ExAC All threshold",
    #        default=100, required=False)
    #parser.add_argument("-XE", "--exacEUR", type=str, dest="exac_eur_threshold", help="ExAC European threshold",
    #        default=100, required=False)
    parser.add_argument("-v", "--verbosity", action="count", default=0)

    args = parser.parse_args()
    out_all = args.out + ".all.bed"
    out_nfe = args.out + ".nfe.bed"

    #outputfile = open(args.out, "w")
    outputfile_all = open(out_all, "w")
    outputfile_nfe = open(out_nfe, "w")

    if args.verbosity >= 2:
        print "{} to the power {} equals {}".format(args.v, args.o, answer)
    elif args.verbosity >= 1:
        print "{}^{} == {}".format(args.x, args.y, answer)

    condel_tbx = pysam.TabixFile("data/fannsdb.small.bed.gz")

    vcf_reader = vcf.Reader(open(args.vcf, 'r'))
    for record in vcf_reader:
        current_chr = "chr" + record.CHROM
        current_id = record.ID
        current_pos = record.POS
        current_ref = record.REF
        #current_filter = ','.join(str(v) for v in record.FILTER)
        current_filter = len(record.FILTER)
        current_alt = ','.join(str(v) for v in record.ALT)

        # SnpEff
        ann = record.INFO['ANN'][0].split('|')
        annotation = ann[1]
        current_gene = ann[3]
        # Condel
        current_condel = getTabixValCondel(condel_tbx, current_chr, current_pos, current_ref, current_alt)

        current_af = ','.join(str(v) for v in record.INFO['AF'])
        ac_adj = ','.join(str(v) for v in record.INFO['AC_Adj'])
        an_adj = record.INFO['AN_Adj']

        #NFE
        ac_nfe = ','.join(str(v) for v in record.INFO['AC_NFE'])
        an_nfe = record.INFO['AN_NFE']

        #AMR
        #nfe_ac = ','.join(str(v) for v in record.INFO['AC_NFE'])
        #nfe_an = record.INFO['AN_NFE']


        out_str = [ current_chr, str(current_pos-1), str(current_pos), current_ref, current_alt, annotation,\
                current_gene, ac_adj , str(an_adj), str(current_condel) ]
        out_nfe = [ current_chr, str(current_pos-1), str(current_pos), current_ref, current_alt, annotation,\
                current_gene, ac_nfe , str(an_nfe), str(current_condel) ]
        out_str = [x or '.' for x in out_str]

        if current_filter == 0:
            #outputfile.write("\t".join(out_str))
            outputfile_all.write("\t".join(out_str))
            outputfile_nfe.write("\t".join(out_nfe))
            outputfile_all.write("\n")
            outputfile_nfe.write("\n")

    #outputfile.close()
    outputfile_all.close()
    outputfile_nfe.close()

if __name__ == "__main__":
    main(sys.argv)

