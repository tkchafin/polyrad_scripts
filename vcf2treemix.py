#!/usr/bin/python

"""
Matt Gibson
Aug. 2019
Indiana University Departement of Biology
Moyle Lab
Converts a vcf file to treemix input. Does not perform any filtering.
Inputs: [1] vcf
        [2] population map
        [3] output name
"""

import sys
import numpy as np
from collections import OrderedDict

def lookup(q, dic):
    rev_dic = {}

    for key, val in dic.items():
        for v in val:
            rev_dic[v] = key
    return(rev_dic[q])


vcf = open(sys.argv[1], 'r')
pops = open(sys.argv[2], 'r')
out1 = open(sys.argv[3], 'w')


popmap = OrderedDict()
for i, line in enumerate(pops):
    l = line.replace('\n','').split()
    popmap[l[0]] = l[1]
    #if l[1] not in popmap.keys():
    #    popmap[l[1]] = [l[0]]
    #else:
    #    popmap[l[1]].append(l[0])
pops.close()

#print(popmap)
#genotypes_bypop = {}
#for pop in popmap.keys():
#    genotypes_bypop[pop] = []



out1.write(' '.join(set(popmap.values())) + '\n')

for i, line in enumerate(vcf):
    l = line.replace("\n", '').split()
    if line[0] == '#' and line[1] == '#':
        next
    elif l[0] == '#CHROM':
        vcf_individuals = l[9:len(l)]
        individual_indexes = [x+9 for x in range(0,len(vcf_individuals))]
        print(vcf_individuals)
        print(individual_indexes)
        assert len(vcf_individuals) == len(individual_indexes)
    else:
        counts = {}
        for k in set(popmap.values()):
            counts[k] = [0,0]
        
        #print(counts)
        
        for i, x in enumerate(individual_indexes):
            call = l[x].split(":")[0]
            #print(call)
            indv = vcf_individuals[i]
            if indv not in popmap.keys():
                continue
            pop = popmap[indv]
            #print(pop)
            if call == "0/0":
                counts[pop][0] += 2
            elif call == "1/1":
                counts[pop][1] += 2
            elif call == "0/1" or call == "1/0":
                counts[pop][0] += 1
                counts[pop][1] += 1
        line_n = []
        #print(counts)
        for k, val in counts.items():
            line_n.append(str(val[0]) + ',' + str(val[1]))
        out1.write(" ".join(line_n) + '\n')
out1.close()
