'''
A tool to calculate two proportion Z-test.
The two input files contain gene variants counts
for their corresponding populations.
Author: Khalid Mahmood
Contact: khalid.mahmood@unimelb.edu.au
Copyright: 2015
'''

#!/usr/bin/python

#from annotations import getTabixValCondel
#from annotations import getExAC
from math import sqrt, erf

import math
import scipy
import scipy.stats
import scipy
import sys
import os
import argparse
import getopt
import vcf
import re
import array
import pysam

'''
def get_score(row, name, action):
    try:
        score = float(row[name])
    except:
        pass
    else:
        action(score)


def read_input(filename):
    gene = []
    acount = []
    anumber = []
    var = []
    individuals = []
    missense_ac	= []
    missense_an = []
    nonsense_ac = []
    nonsense_an = []

    gerp = []
    fathmm = []
    with open(filename) as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            classification = row['benign']
            get_score(row, 'GENE', lambda score: gene.append(score, classification)))
            get_score(row, 'ACount', lambda score: polyphen.append((score, classification)))
            get_score(row, 'ANumber', lambda score: cadd.append((score, classification)))
            get_score(row, 'VAR', lambda score: pri_ph_cons.append((score, classification)))
            get_score(row, 'Individuals', lambda score: gerp.append((score, classification)))
            get_score(row, 'Missense_AC', lambda score: fathmm.append((score, classification)))
            get_score(row, 'Missense_AN', lambda score: fathmm.append((score, classification)))
            get_score(row, 'Nonsense_AC', lambda score: fathmm.append((score, classification)))
            get_score(row, 'Nonsense_AN', lambda score: fathmm.append((score, classification)))

    # sort in ascending order by score
    sift.sort()
    polyphen.sort()
    cadd.sort()
    pri_ph_cons.sort()
    gerp.sort()
    fathmm.sort()

    return sift, polyphen, cadd, pri_ph_cons, gerp, fathmm
'''


def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--A", type=str, dest="file1", help="Input file 1 (tsv)", required=True)
    parser.add_argument("-b", "--B", type=str, dest="file2", help="Input file 2 (tsv)", required=True)
    parser.add_argument("-p", "--P", type=str, dest="p", help="Column number containing the p count", required=True)
    parser.add_argument("-n", "--N", type=str, dest="n", help="Column number containing the sample size n", required=True)
    #parser.add_argument("-X", "--exac", type=str, dest="exac_af_threshold", help="ExAC All threshold",
    #        default=100, required=False)
    #parser.add_argument("-XE", "--exacEUR", type=str, dest="exac_eur_threshold", help="ExAC European threshold",
    #        default=100, required=False)
    parser.add_argument("-v", "--verbosity", action="count", default=0)

    args = parser.parse_args()

    if args.verbosity >= 2:
        print "{} to the power {} equals {}".format(args.v, args.o, answer)
    elif args.verbosity >= 1:
        print "{}^{} == {}".format(args.x, args.y, answer)

    header_str = "Gene\tP1\tT1\tP2\tT2\tPSP\tSE\tz-score\tp-value"

    print header_str

    reader1 = open(args.file1, 'r')
    reader2 = open(args.file2, 'r')
    next(reader1),next(reader2)

    pop1 = {}
    pop2 = {}

    for row in reader1:
        row = row.rstrip('\n')
        row_info = row.split("\t")
        gene = str(row_info[0])
        count = int(row_info[int(args.p)-1])
        number = int(row_info[int(args.n)-1])
        proportion = float(count)/float(number)
        #insert proportion
        #pop1[gene] = round(proportion, 5), number
        pop1[gene] = count, number

    for row in reader2:
        row = row.rstrip('\n')
        row_info = row.split("\t")
        gene = str(row_info[0])
        count = int(row_info[int(args.p)-1])
        number = int(row_info[int(args.n)-1])
        proportion = float(count)/float(number)
        #insert proportion
        #pop2[gene] = round(proportion, 5), number
        pop2[gene] = count, number

    # pooled sample proportion SE Z-score P-value
    for x in pop1:
        if x in pop2:
            count1 = pop1[x][0]
            number1 = pop1[x][1]
            proportion1 = float(count1)/float(number1)

            count2 = pop2[x][0]
            number2 = pop2[x][1]
            proportion2 = float(count2)/float(number2)

            if proportion1>0 and proportion2>0:
                #print str(x) + "\t" + str(pop1[x][0]) + "\t" + str(pop1[x][1]) + "\t" + str(pop2[x][0]) + \
                #    "\t" + str(pop2[x][1]) + "\t",
                print str(x) + "\t" + str(count1) + "\t" + str(number1) + "\t" + str(count2) + \
                    "\t" + str(number2),

                #psp = float(((pop1[x][0]*pop1[x][1])+((pop2[x][0]*pop2[x][1])))/(pop1[x][1]+pop2[x][1]))
                #(D5*157+E5*500)/(157+500)
                psp = float((proportion1*number1)+(proportion2*number2))/float(number1+number2)
                print str(psp) + "\t",

                #standard_error = sqrt(float((psp*(1.0-psp))*(float(1.0/pop1[x][1])+float(1.0/pop2[x][1]))))
                standard_error = sqrt(float(psp*(1.0-psp))*(float(1.0/number1)+float(1.0/number2)))
                print str(standard_error) + "\t",

                #zscore = (pop1[x][0]-pop2[x][0])/standard_error
                zscore = (proportion1-proportion2)/standard_error
                print str(zscore) + "\t",

                #pval = 2*(scipy.stats.norm.sf(abs(zscore)))
                pval = 2*(scipy.stats.norm.sf(abs(zscore)))
                print pval

            #else:
                #print str(x) + "\t" + str(count1) + "\t" + str(number1) + "\t" + str(count2) + \
                #    "\t" + str(number2)
        else:
            print str(x) + " not found."


    reader1.close()
    reader2.close()


if __name__ == "__main__":
    main(sys.argv)

