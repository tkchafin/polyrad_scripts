#!/usr/bin/python

import sys
import os
import getopt

def main():
	params = parseArgs()
	
	if params.popmap:
		popmap = parsePopmap(params.popmap)
	
	#parse VCF
	with open(params.vcf, "r") as vcf:
		last_header=None
		max_ploidy=0
		sample_names=None
		data=list()
		for line in vcf:
			line=line.strip()
			if len(line) == 0:
				continue
			#skip headers
			if line[0] == "#":
				last_header=line.split("\t")
				continue
			else:
				if not sample_names:
					sample_names = last_header[9:]
				fields=line.split("\t")
				
				#check if sample ploidy is higher
				p=max(getPloidies(fields[9:]))
				if p > max_ploidy:
					max_ploidy=p
				
				#add sample data to list
				data.append(fields[9:])
		vcf.close()

	#write output structure
	with open(params.out, "w") as f:
		for i, sample in enumerate(sample_names):
			base= list()
			base.append(sample)
			if params.popmap:
				if sample in popmap.keys():
					base.append(popmap[sample])
				else:
					print("Sample",sample,"not in popmap. Skipping.")
					continue
			if params.extracols:
				for c in range(1, params.extracols):
					base.append("")
			olines=list()
			for a in range(1, max_ploidy):
				olines.append(base)
			for loc in data:
				print(data[i])
				sample_loc = (data[i]).split("/")
				for a in range(1, max_ploidy):
					index=a-1
					if a > len(sample_loc):
						olines[index].append(-9)
					else:
						olines[index].append(sample_loc[index])
				sys.exit()
			for o in olines:
				f.write(o)
				f.write("\n")
		f.close()

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


def getPloidies(samples):
	ploidies=list()
	for s in samples:
		p=len((s.split(":")[0]).split("/"))
		ploidies.append(p)
	return(ploidies)

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hv:o:p:x:', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
		self.out="polyrad.str"
		self.popmap=None
		self.extracols=0
		self.locnames=False


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
			elif opt=="p":
				self.popmap=arg
			elif opt=="x":
				self.extracols=int(arg)
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.vcf:
			self.display_help("Need a VCF file")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\npolyVCFtoStructure.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Converts polyRAD VCF to Structure (.str)")
		print("""
		-v	: VCF input file formatted by polyRAD
		-p	: (Optional) popmap file with pop labels
		-x	: (Optional) number of extra (blank) columns to add
		-o	: Output file name (default=polyrad.vcf)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
