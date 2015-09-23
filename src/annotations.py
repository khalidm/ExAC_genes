'''
A tool to annotate and print variants in tabular format.
Author: Khalid Mahmood (khalid.mahmood@unimelb.edu.au).
Copyright: 2015
'''

#!/usr/bin/python


#from utils import findlist

import sys
import os
import argparse
import getopt
#import vcf
import array
import pysam

#class Error(Exception):
#    """Base-class for exceptions in this module."""

#class UsageError(Error):
#    def __init__(self, msg):
#        self.msg = msg

def getExAC(record, ac_eth, an_eth, index):
    #ac = ','.join(str(v) for v in record.INFO[eth])
    ac = record.INFO[ac_eth][index]
    an = record.INFO[an_eth]
    af = round((float(ac*1.0)/float(an*1.0)),7)
    return str(ac), str(an), str(af)

def getTabixVal(input_tbx, current_chr, current_pos, current_ref, current_alt):
    #current_chr = current_chr.translate(None, 'chr')
    data = input_tbx.fetch(current_chr, current_pos-1, current_pos)
    value = '.'
    if data is not None:
        for row in data:
            row_info = row.split("\t")
            value = row_info[3]
    #else:
    #    value = '.'

    return value

def getTabixValCondel(input_tbx, current_chr, current_pos, current_ref, current_alt):
    #current_chr = current_chr.translate(None, 'chr')
    data = input_tbx.fetch(current_chr, current_pos-1, current_pos)
    value = 0.0001
    if data is not None:
        for row in data:
            row_info = row.split("\t")
            if( current_ref == row_info[3] and current_alt == row_info[4] ):
                value = row_info[7]
                break
    return round(float(value), 4)

def getfathmm(fathmm_tbx, current_chr, current_pos, current_ref, current_alt):
    #current_chr = current_chr.translate(None, 'chr')
    data = fathmm_tbx.fetch(current_chr, current_pos-1, current_pos)
    fathmm_score = 0.0
    if data is not None:
        for row in data:
            row_info = row.split("\t")
            fathmm_ref = row_info[3]
            fathmm_alt = row_info[4]
            if(fathmm_ref == current_ref and fathmm_alt == current_alt):
                fathmm_score = row_info[5]
                break

    return fathmm_score

def getTabixBool(input_tbx, current_chr, current_pos, current_ref, current_alt):
    #current_chr = current_chr.translate(None, 'chr')
    data = input_tbx.fetch(current_chr, current_pos-1, current_pos)
    val = '.'

    if data is not None:
        for row in data:
            #print current_chr + ":" + str(current_pos) + ":" + str(row.split("\t"))
            val = "T"

    return val

def adjust_scores(condel, sift, polyphen, fathmm, annotation):

    if( ("stop_lost" in annotation) or ("stop_gained" in annotation) ):
        condel = 0.9999
        sift = 0.0001
        polyphen = 0.9999
        fathmm = 0.9999

    return condel, sift, polyphen, fathmm

