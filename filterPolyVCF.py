#!/usr/bin/python

import sys
import os
import getopt
import random

def main():
	params = parseArgs()
	
	f = open(params.out, 'w')
	
	with open(params.vcf, "r") as vcf:
		this_loc=None
		loc_snps=list()
		for line in vcf:
			line=line.strip()
			#directly transfer header lines
			if line[0] == "#":
				f.write(line)
				f.write("\n")
			else:
				fields=line.split("\t")
				if not this_loc:
					this_loc=fields[0]
				else:
					if fields[0] != this_loc:
						#submit locus here
						if params.random:
							#randomly select a SNP to write
							if len(this_loc)>1:
								f.write("\t".join(random.choice(this_loc)))
								f.write("\n")
							else:
								f.write("\t".join(this_loc[0]))
								f.write("\n")
						else:
							for s in loc_snps:
								f.write("\t".join(s))
								f.write("\n")
						#clear locus data
						this_loc=fields[0]
						loc_snps=list()
						
				ref=fields[3].split(",")
				alt=fields[4].split(",")
				
				#if biallelic filter on and site has >2 alleles, skip
				if params.biallelic == True:
					if (len(ref) + len(alt) > 2):
						continue
		
				#grab info fields if field filters are turned on
				if params.noFilter == False:
					info=parseINFO(fields[7])
				
					#filter on HH
					if params.minH:
						if "HH" in info.keys():
							if info["HH"] < params.minH:
								continue #skip locus
						else:
							print("ERROR: Key \"HH\" not in INFO field")
							print(info)
							sys.exit()
					if params.maxH:
						if "HH" in info.keys():
							if info["HH"] > params.maxH:
								continue #skip locus
						else:
							print("ERROR: Key \"HH\" not in INFO field")
							print(info)
							sys.exit()
					if params.minC:
						if "NS" in info.keys():
							if info["NS"] < params.minC:
								continue #skip locus
						else:
							print("ERROR: Key \"NS\" not in INFO field")
							sys.exit()
					if params.minD:
						if "DP" in info.keys():
							if info["DP"] < params.minD:
								continue #skip locus
						else:
							print("ERROR: Key \"DP\" not in INFO field")
							sys.exit()
				if params.minP or params.maxP:
					p=getPloidies(fields[9:])
					minP=min(p)
					maxP=max(p)
					if params.minP and minP < params.minP:
						continue
					if params.maxP and maxP > params.maxP:
						continue
				
				#if all passed, append to locus list
				loc_snps.append(fields)

		vcf.close()
	f.close()

def getPloidies(samples):
	ploidies=list()
	for s in samples:
		p=len((s.split(":")[0]).split("/"))
		ploidies.append(p)
	return(ploidies)
		

def parseINFO(info):
	d=dict()
	l=info.split(";")
	for f in l:
		fl=f.split("=")
		d[fl[0]]=float(fl[1])
	return(d)

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hv:o:bc:p:P:m:M:rd:', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
		self.out="polyrad.vcf"
		self.biallelic=False
		self.minH=None
		self.maxH=None
		self.minC=None
		self.minP=None
		self.maxP=None
		self.random=False
		self.noFilter=True
		self.minD=None

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
			elif opt=="v":
				self.vcf=arg
			elif opt=="o":
				self.out=arg
			elif opt=="b":
				self.biallelic=True
			elif opt=="p":
				if int(arg)<1:
					self.display_help("-p cannot be less than 1")
				self.minP=int(arg)
			elif opt=="P":
				self.maxP=int(arg)
			elif opt=="m":
				if float(arg) < 0.0:
					self.display_help("-m cannot be less than 1")
				self.noFilter=False
				self.minH=float(arg)
			elif opt=="M":
				self.noFilter=False
				self.maxH=float(arg)
			elif opt=="c":
				self.noFilter=False
				self.minC=int(arg)
			elif opt=="d":
				self.minD=int(arg)
				self.noFilter=False
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.vcf:
			self.display_help("Need an ipyrad VCF file")
		if self.maxP and self.minP > self.maxP:
			self.display_help("-p can't be greater than -P")
		if self.maxH < self.minH:
			self.display_help("-m can't be greater than -M")
		
		if self.noFilter:
			if not self.biallelic and not self.random:
				if not self.minP and not self.maxP:
					self.display_help("No filters set. Nothing to do.")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nfilterPolyVCF.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description:Filters the VCF file produced by polyRAD's RADdata2VCF() function")
		print("""
	Arguments:
		-v	: VCF file
		-b	: (Boolean) Toggle to only retain [default=off]
		-c	: Mimumum number of samples genotyped (NS) [default=None]
		-p	: Minimum allowable ploidy per locus [default=None]
		-P	: Maximum allowable ploidy per locus [default=None]
		-m	: Minimum HindHe/ HH value [default=None]
		-M	: Maximum HindHe/ HH value [default=None]
		-r	: (Boolean) Randomly sample 1 SNP per locus [default=off]
		-d	: Minimum DP [default=None]
		-o	: Output file name (default=polyrad.vcf)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
