#!/usr/bin/python

import sys
import os
import getopt

def main():
	params = parseArgs()
	
	f = open(params.out, 'w')
	
	with open(params.vcf, "r") as vcf:
		last_header=None
		for line in vcf:
			line=line.strip()
			#skip headers
			if line[0] == "#":
				last_header=line.split("\t")
				continue
			else:
				fields=line.split("\t")

		vcf.close()
	f.close()


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
\		-o	: Output file name (default=polyrad.vcf)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
