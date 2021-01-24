import sys
import os
import numpy as np

# """
# Implements the lnRH statistics of Schlötterer & Dieringer (2005)
# A Novel Test Statistic for the Identification of Local Selective Sweeps Based on Microsatellite Gene Diversity
# """
# def lnRH(pi1, pi2):
# 	rh=0.0
# 	he1=1.0-(polyHe(pi1))
# 	he2=1.0-(polyHe(pi2))
# 
# 	if he2 == 1.0 or he1 == 1.0:
# 		return(0.0)
# 	rh=np.log((np.square(1.0/(1.0-he1))-1)/(np.square(1.0/(1.0-he2))-1))
# 	return(rh)

"""
Calculates locus Dxy as:
Pi(x)*(1-Pi(y)) + Pi(y)*(1-Pi(x))
Warning: Assumes locus is bi-allelic!
"""
def Dxy(pi1, pi2):
	Dxy=0.0
	pix=0.0
	piy=0.0
	if "0" in pi1.keys():
		pix=pi1['0']
	else:
		pix=0.0
	if "0" in pi2.keys():
		piy=pi2['0']
	else:
		piy=0.0
	Dxy=float((pix*(1.0-piy))+(piy*(1.0-pix)))
	# print("piy:",piy)
	# print("pix:",pix)
	# print("Dxy:", Dxy)
	return(Dxy)
		
"""
Computes the D statistic for Jost (2008)
GST and its relatives do not measure differentiation
expressed in terms of He (or Hs) and Ht 
See Meirmans et al (2018)
"""
def JostD(pi1, pi2):
	dest=0.0
	he1=polyHe(pi1)
	he2=polyHe(pi2)
	Hs=(he1+he2)/2.0
	Ht=getHt(pi1, pi2)
	#print("Hs: ", Hs)
	#print("Ht: ", Ht)
	dest=2.0*(float(Ht-Hs)/(1.0-Hs))
	#print("D: ", dest)
	return(dest)

"""
G''ST following Meirmans et al 2018
"""
def GppST(pi1, pi2):
	he1=polyHe(pi1)
	he2=polyHe(pi2)
	Hs=(he1+he2)/2.0
	Ht=getHt(pi1, pi2)
	if Ht==Hs:
		return(0.0)
	num=2.0*(Ht-Hs)
	den=(2.0*Ht-Hs)*(1.0-Hs)
	# print("Hs:",Hs)
	# print("Ht:",Ht)
	# print("num:",num)
	# print("den:",den)
	return(num/den)

def getHt(pi1, pi2):
	Ht=0.0
	alleles=list(pi1.keys())
	alleles.extend(list(pi2.keys()))
	for allele in set(alleles):
		f1=0.0
		f2=0.0
		if allele in pi1.keys():
			f1=pi1[allele]
		if allele in pi2.keys():
			f2=pi2[allele]
		m=float(f1+f2)/2.0
		#print(m)
		Ht+=(m*m)
		#print(Ht)
	#print(Ht)
	return(1.0-Ht)

def isMonomorphic(pops):
	for pop in pops.keys():
		p=countsFromGenotypes(pops[pop], Pi=True, skipMissing=True)
		if len(p.keys()) > 1:
			return(False)
	return(True)

def isBiallelic(pops):
	all=list()
	for pop in pops.keys():
		all.extend(pops[pop])
	p=countsFromGenotypes(all, Pi=True, skipMissing=True)
	if len(p.keys())==2:
		return(True)
	else:
		return(False)

"""
Calculates a gene diversity measure acceptable for polyploids 
Same as expected heterozygosity in diploids 
Following recommendation of Meirmans et al (2018)
The Analysis of Polyploid Genetic Data
"""
def polyHe(pi):
	he=0.0
	sumPiSq=0.0
	for nuc in pi.keys():
		sumPiSq+=float(pi[nuc]*pi[nuc])
	he=(sumPiSq)
	return(1.0-he)


def countsFromGenotypes(inds, Pi=False, skipMissing=True):
	c=dict()
	decomp=list()
	total=0.0
	for i in inds:
		decomp.extend(i.split("/"))
	for nuc in set(decomp):
		if nuc in ["-", "?", "N", "."]:
			if not skipMissing:
				total+=decomp.count(nuc)
			continue
		else:
			c[nuc]=decomp.count(nuc)
			total+=decomp.count(nuc)
	if Pi:
		for nuc in c.keys():
			c[nuc] = float(c[nuc]/total)
	return(c)