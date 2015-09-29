'''
A tool extract gene/burden from ExAC vcf files
Author: Khalid Mahmood
Contact: khalid.mahmood@unimelb.edu.au
Copyright: 2015
'''

#!/usr/bin/python

from annotations import getTabixValCondel
from annotations import getExAC

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
    out_combined = args.out + ".bed"
    out_all = args.out + ".all.bed"
    out_nfe = args.out + ".nfe.bed"
    out_fin = args.out + ".fin.bed"
    out_afr = args.out + ".afr.bed"
    out_eas = args.out + ".eas.bed"
    out_sas = args.out + ".sas.bed"
    out_amr = args.out + ".amr.bed"
    out_oth = args.out + ".oth.bed"

    outputfile = open(out_combined, "w")
    outputfile_all = open(out_all, "w")
    outputfile_nfe = open(out_nfe, "w")
    outputfile_fin = open(out_fin, "w")
    outputfile_afr = open(out_afr, "w")
    outputfile_eas = open(out_eas, "w")
    outputfile_sas = open(out_sas, "w")
    outputfile_amr = open(out_amr, "w")
    outputfile_oth = open(out_oth, "w")

    if args.verbosity >= 2:
        print "{} to the power {} equals {}".format(args.v, args.o, answer)
    elif args.verbosity >= 1:
        print "{}^{} == {}".format(args.x, args.y, answer)

    condel_tbx = pysam.TabixFile("data/fannsdb.small.bed.gz")

    header_str = "#CHR\tSTART\tEND\tREF\tALT\tANNOTATION\tGENE\tAC\tAN\tAF\tSAMPLES\tCONDEL\n"

    outputfile_all.write(header_str),outputfile_nfe.write(header_str)
    outputfile_fin.write(header_str),outputfile_afr.write(header_str)
    outputfile_eas.write(header_str),outputfile_sas.write(header_str)
    outputfile_amr.write(header_str),outputfile_oth.write(header_str)

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
        #ann = record.INFO['ANN'].split(',')
        #ann = record.INFO['ANN'][0].split('|')
        #annotation = ann[1]
        #current_gene = ann[3]
        # Condel
        condels = []
        alt_index = 0

        if current_filter == 0:

            # Iteratre though each ALT allele
            for v in record.ALT:
                if len(record.INFO['ANN']) < len(record.ALT):
                    ann_temp = record.INFO['ANN'][0].split('|')
                    annotation = ann_temp[1]
                    current_gene = ann_temp[3]
                else:
                    ann_temp = record.INFO['ANN'][alt_index].split('|')
                    if v==ann_temp[0]:
                        annotation = ann_temp[1]
                        current_gene = ann_temp[3]
                    else:
                        annotation = ann_temp[1]
                        current_gene = ann_temp[3]

                #print str(current_chr) + "\t" + str(current_pos) + "\t" + str(current_ref) + "\t" + str(v)
                #condels.append(str(getTabixValCondel(condel_tbx, current_chr, current_pos, current_ref, v)))

                # current_condel = ','.join(str(p) for p in condels)
                # current_af = ','.join(str(v) for v in record.INFO['AF'])

                if any(x in ['frameshift_variant', 'start_lost', 'stop_gained'] for x in annotation.split('&')):
                    current_condel = 1.000
                else:
                    current_condel = getTabixValCondel(condel_tbx, current_chr, current_pos, current_ref, v)

                # print current_chr + "\t" + str(current_pos) + " fileter = " + str(current_filter)

                ac_adj, an_adj, af_adj, all_samples = getExAC(record, "AC_Adj", "AN_Adj", alt_index)
                ac_nfe, an_nfe, af_nfe, nfe_samples = getExAC(record, "AC_NFE", "AN_NFE", alt_index)
                ac_fin, an_fin, af_fin, fin_samples = getExAC(record, "AC_FIN", "AN_FIN", alt_index)
                ac_afr, an_afr, af_afr, afr_samples = getExAC(record, "AC_AFR", "AN_AFR", alt_index)
                ac_eas, an_eas, af_eas, eas_samples = getExAC(record, "AC_EAS", "AN_EAS", alt_index)
                ac_sas, an_sas, af_sas, sas_samples = getExAC(record, "AC_SAS", "AN_SAS", alt_index)
                ac_amr, an_amr, af_amr, amr_samples = getExAC(record, "AC_AMR", "AN_AMR", alt_index)
                ac_oth, an_oth, af_oth, oth_samples = getExAC(record, "AC_OTH", "AN_OTH", alt_index)

                alt_index = alt_index + 1

                out_combined_str = [ current_chr, str(current_pos-0), str(current_pos), current_ref, str(v), current_gene,\
                        annotation, str(current_condel), af_adj, af_nfe, af_fin, af_afr, af_eas, af_sas, af_amr, af_oth ]
                out_all = [ current_chr, str(current_pos-0), str(current_pos), current_ref, str(v), annotation,\
                        current_gene, ac_adj , an_adj, af_adj, all_samples, str(current_condel) ]
                out_nfe = [ current_chr, str(current_pos-0), str(current_pos), current_ref, str(v), annotation,\
                        current_gene, ac_nfe , an_nfe, af_nfe, nfe_samples, str(current_condel) ]
                out_fin = [ current_chr, str(current_pos-0), str(current_pos), current_ref, str(v), annotation,\
                        current_gene, ac_fin, an_fin, af_fin, fin_samples, str(current_condel) ]
                out_afr = [ current_chr, str(current_pos-0), str(current_pos), current_ref, str(v), annotation,\
                        current_gene, ac_afr , an_afr, af_afr, afr_samples, str(current_condel) ]
                out_eas = [ current_chr, str(current_pos-0), str(current_pos), current_ref, str(v), annotation,\
                        current_gene, ac_eas , an_eas, af_eas, eas_samples, str(current_condel) ]
                out_sas = [ current_chr, str(current_pos-0), str(current_pos), current_ref, str(v), annotation,\
                        current_gene, ac_sas , an_sas, af_sas, sas_samples, str(current_condel) ]
                out_amr = [ current_chr, str(current_pos-0), str(current_pos), current_ref, str(v), annotation,\
                        current_gene, ac_amr , an_amr, af_amr, amr_samples, str(current_condel) ]
                out_oth = [ current_chr, str(current_pos-0), str(current_pos), current_ref, str(v), annotation,\
                        current_gene, ac_oth, an_oth, af_oth, oth_samples, str(current_condel) ]

                out_combined_str = [x or '.' for x in out_combined_str]
                out_all = [x or '.' for x in out_all]
                out_nfe = [x or '.' for x in out_nfe]
                out_fin = [x or '.' for x in out_fin]
                out_afr = [x or '.' for x in out_afr]
                out_eas = [x or '.' for x in out_eas]
                out_sas = [x or '.' for x in out_sas]
                out_amr = [x or '.' for x in out_amr]
                out_oth = [x or '.' for x in out_oth]

                outputfile.write("\t".join(out_combined_str)),outputfile.write("\n")
                outputfile_all.write("\t".join(out_all)),outputfile_nfe.write("\t".join(out_nfe))
                outputfile_fin.write("\t".join(out_fin)),outputfile_afr.write("\t".join(out_afr))
                outputfile_eas.write("\t".join(out_fin)),outputfile_sas.write("\t".join(out_afr))
                outputfile_amr.write("\t".join(out_fin)),outputfile_oth.write("\t".join(out_afr))
                outputfile_all.write("\n"),outputfile_nfe.write("\n")
                outputfile_fin.write("\n"),outputfile_afr.write("\n")
                outputfile_eas.write("\n"),outputfile_sas.write("\n")
                outputfile_amr.write("\n"),outputfile_oth.write("\n")

    outputfile.close()
    outputfile_all.close(),outputfile_nfe.close()
    outputfile_fin.close(),outputfile_afr.close()
    outputfile_eas.close(),outputfile_sas.close()
    outputfile_amr.close(),outputfile_oth.close()

if __name__ == "__main__":
    main(sys.argv)

