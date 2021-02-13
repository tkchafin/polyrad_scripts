#!/usr/bin/python

import sys
import os
import getopt
import pandas as pd
import operator
import polyStats as ps

def main():
	params = parseArgs()
	
	popmap=parsePopmap(params.popmap)
	
	popdicts=getPopDicts(popmap, params.pop1, params.pop2)
	
	name=list()
	pos=list()
	fA=list()
	fB=list()
	
	with open(params.vcf, "r") as vcf:
		this_loc=None
		sample_indices=dict()
		offset=9
		first=True
		last_line=list()
		for line in vcf:
			line=line.strip()
			#directly transfer header lines
			if line[0] == "#":
				last_line=line.split("\t")
			else:
				#if first data line, grab sample names from last header line
				if first:
					sample_names=last_line[offset:]
					index=0
					for samp in sample_names:
						sample_indices[samp]=offset+index
						index+=1
					first=False

				fields=line.split("\t")
				locus=fields[0]
				position=fields[1]
				
				gen = getGenotypes(fields, popdicts, sample_indices)
				
				#print(gen)
				
				pi=getPis(gen)
				
				maxA = max(pi["pop1"].items(), key=operator.itemgetter(1))[0]
				
				#print(maxA)
				afreq=pi["pop1"][maxA]
				bfreq=0.0
				if maxA in pi["pop2"].keys():
					bfreq=pi["pop2"][maxA]
				
				fA.append(afreq)
				fB.append(bfreq)

				
				name.append(locus)
				pos.append(position)

		vcf.close()
	stuff = {'Locus': name, 'Position': pos, 'FreqA': fA,
	'FreqB': fB}
	df = pd.DataFrame(stuff)
	print(df)
	df.to_csv(params.oname, sep="\t", header=True, index=False)  

def getPis(gen):
	pis=dict()
	for pop in gen.keys():
		if pop not in pis.keys():
			pis[pop]=list()
		pis[pop]= ps.countsFromGenotypes(gen[pop], Pi=True, skipMissing=True)
	genT=gen["pop1"]
	genT.extend(gen["pop2"])
	piT=ps.countsFromGenotypes(genT, Pi=True, skipMissing=True)
	pis["popT"]=piT #Ht
	return(pis)

def getGenotypes(fields, popdicts, sample_indices):
	gen=dict()
	for pop in popdicts.keys():
		if pop not in gen.keys():
			gen[pop]=list()
		for ind in popdicts[pop]:
			if ind in sample_indices.keys():
				gt=fields[sample_indices[ind]]
				gen[pop].append(gt.split(":")[0])
	return(gen)
	
	
def getPopDicts(popmap, pop1, pop2):
	popdicts=dict()
	pop2d=make2Dpopmap(popmap)
	
	popdicts["pop1"]=list()
	popdicts["pop2"]=list()
	for pop in pop1:
		popdicts["pop1"].extend(pop2d[pop])
	for pop in pop2:
		popdicts["pop2"].extend(pop2d[pop])

	return(popdicts)

#Makes a dict of lists from a popmap
def make2Dpopmap(p):
	ret = dict()
	for s in p:
		if p[s] not in ret:
			ret[p[s]] = list()
		ret[p[s]].append(s)
	return(ret)

#function reads a tab-delimited popmap file and return dictionary of assignments
def parsePopmap(popmap):

	ret = dict()
	if os.path.exists(popmap):
		with open(popmap, 'r') as fh:
			try:
				contig = ""
				seq = ""
				for line in fh:
					line = line.strip()
					if not line:
						continue
					else:
						stuff = line.split()
						ret[stuff[0]] = stuff[1]
				return(ret)
			except IOError:
				print("Could not read file ",pairs)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%popmap)

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'h1:2:p:v:o:', \
			["help", "popmap=", "vcf=", "pop1=", "pop2=", "oname="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
		self.popmap=None
		self.pop1=None
		self.pop2=None
		self.oname="out.tsv"


		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == "h" or opt == "help":
				continue
			elif opt=="vcf" or opt=="v":
				self.vcf=arg
			elif opt=="pop1" or opt=="1":
				self.pop1=arg.split("+")
			elif opt=="pop2" or opt=="2":
				self.pop2=arg.split("+")
			elif opt=="popmap" or opt=="p":
				self.popmap=arg
			elif opt=="o" or opt=="oname":
				self.oname=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.vcf:
			self.display_help("No VCF file provided.")
		if not self.popmap:
			self.display_help("No popmap file provided.")
		if not self.pop1 or not self.pop2:
			self.display_help("Populations not specified")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nfreqAB.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Calculates frequency of Population 1 major allele in Population 2")
		print("""
		-v,--vcf	: Input polyVCF file 
		-p,--popmap	: Tab-delimitation population map file
		-1,--pop1	: Identifier for population 1 
			NOTE: multiple can be separated by + (e.g. popA+popB)
		-2,--pop2	: Identifier for population 2
		-o,--oname	: Output file name [default=out.tsv]
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
