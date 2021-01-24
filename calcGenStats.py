#!/usr/bin/python

import sys
import os
import getopt
import pandas as pd

import polyStats as ps

def main():
	params = parseArgs()
	
	popmap=parsePopmap(params.popmap)
	
	popdicts=getPopDicts(popmap, params.pop1, params.pop2, params.pop3)
	
	name=list()
	pos=list()
	dxy=list()
	dxo=list()
	dyo=list()
	rnd=list()
	gpst=list()
	jost=list()
	isMon=list()
	isBi=list()
	
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
				
				mon=ps.isMonomorphic(gen)
				#print("Monomorphic?:",mon)
				
				bi=ps.isBiallelic(gen)
				#print("Biallelic?:",bi)
				
				# rh=ps.lnRH(pi["pop1"], pi["pop2"])
				# print("lnRH:",rh)
				
				#calculate Dxy, Dxo, Dxy, RND(Dxy/mean(Dxo,Dxy))
				Dxy=ps.Dxy(pi["pop1"], pi["pop2"])
				#print("Dxy:",Dxy)
				
				Dxo=ps.Dxy(pi["pop1"], pi["out"])
				#print("Dxo:",Dxo)
				
				Dyo=ps.Dxy(pi["pop2"], pi["out"])
				#print("Dyo:",Dyo)
				
				if Dxo == 0.0 and Dyo==0.0:
					RND=0.0
				else:
					RND=Dxy/((Dxo+Dyo)/2.0)
				#print("RND:",RND)
				
				#G''st (Meirmans and Hedrick 2011)
				Gst=ps.GppST(pi["pop1"], pi["pop2"])
				#print("G''st:", Gst)
				
				#Dest (Jost 2008)
				D=ps.JostD(pi["pop1"], pi["pop2"])
				#print("D:",D)
				
				name.append(locus)
				pos.append(position)
				gpst.append(Gst)
				jost.append(D)
				isMon.append(mon)
				isBi.append(bi)
				dxy.append(Dxy)
				dxo.append(Dxo)
				dyo.append(Dyo)
				rnd.append(RND)
		vcf.close()
	stuff = {'Locus': name, 'Position': pos, 'Gst': gpst,
	'JostD': jost, 'Dxy': dxy, 'Dxo': dxo, 'Dyo': dyo,
	'RND': rnd, 'Monomorphic':isMon, 'Biallelic':isBi}
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
	piT=ps.countsFromGenotypes(genT)
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
	
	
def getPopDicts(popmap, pop1, pop2, pop3):
	popdicts=dict()
	pop2d=make2Dpopmap(popmap)
	
	popdicts["pop1"]=list()
	popdicts["pop2"]=list()
	popdicts["out"]=list()
	for pop in pop1:
		popdicts["pop1"].extend(pop2d[pop])
	for pop in pop2:
		popdicts["pop2"].extend(pop2d[pop])
	
	if pop3 is None:
		#ge all pops not in pop1 or pop2
		for pop in pop2d.keys():
			if pop not in pop1 and pop not in pop2:
				popdicts["out"].extend(pop2d[pop])
	else:
		for pop in pop3:
			popdicts["out"].extend(pop2d[pop])

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
			options, remainder = getopt.getopt(sys.argv[1:], 'h1:2:3:p:v:o:', \
			["help", "popmap=", "vcf=", "pop1=", "pop2=", "pop3=", "oname="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
		self.popmap=None
		self.pop1=None
		self.pop2=None
		self.pop3=None
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
			elif opt=="pop3" or opt=="3":
				self.pop3=arg.split("+")
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
		print ("\ncalGenStats.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Calculated RNDmin, Dxy, and lnRH from a polyVCF")
		print("""
		-v,--vcf	: Input polyVCF file 
		-p,--popmap	: Tab-delimitation population map file
		-1,--pop1	: Identifier for population 1 
			NOTE: multiple can be separated by + (e.g. popA+popB)
		-2,--pop2	: Identifier for population 2
		-3,--pop3	: (Optional) outgroup for RND
		-o,--oname	: Output file name [default=out.tsv]
			If no outgroups provided, script will use all other samples for RNDmin
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
