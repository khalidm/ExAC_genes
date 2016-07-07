import sys
import os
import argparse

from string import letters
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats

pd.options.mode.chained_assignment = None

def sort_props(bd_types, filter_p1):
    #sdfs
    current_file = ''
    types_length = len(bd_types)

    #bd2 = [ '0.5+p', '0.6+p', '0.7+p', '0.5c', '0.6c', '0.7c' ]
    #output_matrix = np.zeros(shape=(types_length, types_length))
    #output_matrix = pd.DataFrame(index=bd2,columns=bd2)
    #output_matrix = output_matrix.fillna(0.0)

    gene_lists = []

    for i in bd_types:
        current_file = "proportions" + "/" + i + ".prop.txt"
        data = pd.read_csv(current_file, sep='\t', header=0)
        data = data[data['pvalue']<0.01]
        data = data[data['p1']>filter_p1]
        #data = data[data['pvalue']<0.001]
        data = data.sort_values(['p1'], ascending=0)
        tp = type(data)
        gene_lists.append(data.gene.tolist())
        outputfilename = "proportions" + "/" + i + ".prop.filtered.txt"
        data.to_csv(outputfilename, sep='\t')

    

def get_ktau_scores(bd_types, top, filter_p1):
    current_file = ''
    types_length = len(bd_types)

    bd2 = [ '0.5CP', '0.6CP', '0.7CP', '0.5C', '0.6C', '0.7C' ]
    
    #output_matrix = np.zeros(shape=(types_length, types_length))
    output_matrix = pd.DataFrame(index=bd2,columns=bd2)
    output_matrix = output_matrix.fillna(0.0)


    gene_lists = []

    for i in bd_types:
        current_file = "proportions" + "/" + i + ".prop.txt"
        data = pd.read_csv(current_file, sep='\t', header=0)
        data = data[data['pvalue']<0.1]
        data = data[data['p1']>filter_p1]
        #data = data[data['pvalue']<0.001]
        data = data.sort_values(['p1'], ascending=0)
        tp = type(data)
        gene_lists.append(data.gene.tolist())

    current_file = ''

    i = 0
    j = 0
    #output_matrix.ix['0.5+p']['0.6+p'] = 1
    #print gene_lists[0].shape()
    #print gene_lists[1].shape()

    for i in range(types_length):
        for j in range(types_length):
            top1 = len(gene_lists[i])
            top2 = len(gene_lists[j])
            top = min(top1, top2)/2
            tau, pval = stats.kendalltau(gene_lists[i][0:top], gene_lists[j][0:top])
            #tau, pval = stats.spearmanr(gene_lists[i][0:top], gene_lists[j][0:top])
            #tau, pval = stats.ranksums(gene_lists[i][0:top], gene_lists[j][0:top])
            #print str(tau) + ":" + str(pval)
            if(pval<10):
                output_matrix.ix[i][j] = tau

    return output_matrix

def get_intersection(bd_types, filter_p1):
    current_file = ''
    types_length = len(bd_types)

    bd2 = [ '0.5CP', '0.6CP', '0.7CP', '0.5C', '0.6C', '0.7C' ]
    
    # intersection
    output_intersection_matrix = pd.DataFrame(index=bd2,columns=bd2)
    output_intersection_matrix = output_intersection_matrix.fillna(0.0)

    gene_lists = []

    for i in bd_types:
        current_file = "proportions" + "/" + i + ".prop.txt"
        data = pd.read_csv(current_file, sep='\t', header=0)
        data = data[data['pvalue']<0.1]
        data = data[data['p1']>filter_p1]
        #data = data[data['pvalue']<0.001]
        data = data.sort_values(['p1'], ascending=0)
        tp = type(data)
        gene_lists.append(data.gene.tolist())

    current_file = ''

    i = 0
    j = 0
    #output_matrix.ix['0.5+p']['0.6+p'] = 1
    #print gene_lists[0].shape()
    #print gene_lists[1].shape()

    for i in range(types_length):
        for j in range(types_length):
            top1 = len(gene_lists[i])
            top2 = len(gene_lists[j])
            top = min(top1, top2)
            x = len(np.intersect1d(gene_lists[i][0:top], gene_lists[j][0:top]))
            int_proportion = x/float(top)
            #print str(int_proportion)
            if(x>=0):
                output_intersection_matrix.ix[i][j] = float(int_proportion)
                #output_intersection_matrix.ix[i][j] = int(x)
    
    return output_intersection_matrix 


def main(argv):
    #sdfsdf
    bd_types = [ 'bd10', 'bd20', 'bd30', 'bd10C', 'bd20C', 'bd30C', 'bd100', 'bd200', 'bd300', 'bd100C', 'bd200C', 'bd300C' ]
    bd_types0 = [ 'bd10', 'bd20', 'bd30', 'bd10C', 'bd20C', 'bd30C' ] # , 'bd100', 'bd200', 'bd300', 'bd100C', 'bd200C', 'bd300C' ]
    bd_types00 = [ 'bd100', 'bd200', 'bd300', 'bd100C', 'bd200C', 'bd300C' ]

    top=500
    filter_p1 = 0

    sort_props(bd_types, filter_p1)

    sort_mat_output0 = sort_props(bd_types0, filter_p1)
    sort_mat_output00 = sort_props(bd_types00, filter_p1)

    out_mat_0 = get_ktau_scores(bd_types0, top, filter_p1)
    out_mat_00 = get_ktau_scores(bd_types00, top, filter_p1)

    output_intersection_matrix0 = get_intersection(bd_types0, filter_p1)
    output_intersection_matrix00 = get_intersection(bd_types00, filter_p1)
    

    sns.set(style="white")
    sns.set(font="monospace", font_scale=0.4)

    #mask = np.zeros_like(out_mat_0, dtype=np.bool)
    #mask[np.triu_indices_from(mask)] = True
    mask0 = np.tri(out_mat_0.shape[0], k=-1)
    mask00 = np.tri(out_mat_00.shape[0], k=-1)

    #print out_mat_0
    fig, ax = plt.subplots(2,2)
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    ax[0][0] = sns.heatmap(out_mat_0, mask = mask0, cmap=cmap, vmax=.3, square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax = ax[0][0], annot=True)
    ax[0][1] = sns.heatmap(out_mat_00, mask = mask00, cmap=cmap, vmax=.3, square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax = ax[0][1], annot=True)
    sns.set(font="monospace", font_scale=0.25)
    ax[1][0] = sns.heatmap(output_intersection_matrix0, mask = mask0, square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax = ax[1][0], annot=True)
    ax[1][1] = sns.heatmap(output_intersection_matrix00, mask = mask00, square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax = ax[1][1], annot=True)
    
    #sns.heatmap(out_mat_0, mask=mask, cmap=cmap, vmax=.3, square=True, xticklabels=5, yticklabels=5, linewidths=.5, cbar_kws={"shrink": .5}, ax = ax)
    #ax.set(xlabel='common xlabel', ylabel='common ylabel')    
    #fig = ax.get_figure()
    fig.savefig("output.pdf")


    # calculate kendall's tau
    #index = 0
    #for x in bd_types:
    #    tau, p_value = stats.kendalltau(gene_lists[index][0:500], gene_lists[0][0:5])
    #    index = index + 1

    #print gene_lists[0][0:5]
    #print gene_lists[1][0:5]
    #x1 = [12, 2, 1, 12, 2]
    #x2 = [1, 4, 7, 1, 0]

    #tau, p_value = stats.kendalltau(x1, x2)
    #print str(tau) + ":" + str(p_value)


if __name__ == "__main__":
    main(sys.argv)

