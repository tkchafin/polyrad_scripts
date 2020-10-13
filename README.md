# polyrad_scripts
Scripts for working with inputs and outputs of [polyRAD](https://github.com/lvclark/polyRAD) genotyper. If you find these useful please cite the paper(s) describing polyRAD: [Clark et al 2019](https://doi.org/10.1534/g3.118.200913) and [Clark et al 2020](https://doi.org/10.1101/2020.01.11.902890).

## ipyrad2polyrad.py

Python script for converting the VCF file output by [ipyrad](https://ipyrad.readthedocs.io/en/latest/) to the proper format for importing into polyRAD. This is necessary because the allele depths are not formatted properly as output by ipyrad. 

```
./ipyrad2polyrad.py -h
Exiting because help menu was called.

ipyrad2polyrad.py

Author: Tyler K Chafin, University of Arkansas
Contact: tkchafin@uark.edu
Description:Converts the ipyrad VCF to a format usable for polyRAD

		-v,--vcf	: VCF input with ipyrad "CATG" field
		-b,--biallelic	: [Boolean] Toggle to skip non-biallelic sites
		-o,--out	: Output file name (default=polyrad.vcf)
```
NOTE: If using ipyrad for processing RAD data prior to genotyping in polyRAD, be sure to relax the parameters in ipyrad aimed at removing merged paralogs (i.e. max shared heterozygosity and max alleles per individual), otherwise these will not be present in the output VCF file!

## filterHHvcf.py

Filtering the VCF file output by polyRAD's RADdata2VCF() function using the Hind/He statistic (make sure to set hinde=TRUE in RADdata2VCF)

```
./filterHHvcf.py -h
Coming soon
```

## polyVCFtoStructure.py 

This script converts a VCF file containing polyploid or mixed ploidy genotypes to a specially formatted file for analysis of population structure in the program [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html). 

In this file, genotypes will be formatted as X lines per individual, where X is the maximum ploidy detected in the input VCF file. Any individuals (or genotypes, if ploidy varies by locus within individuals) will be recorded as 'missing data' (-9) for unneeded fields, following [Stift et al 2019](https://www.nature.com/articles/s41437-019-0247-6).

For example, in the following file IndivA is diploid, Indiv2 is tetraploid, and Indiv3 has mixed-ploidy genotypes:
```
Indiv1  0 1 1 0 -9
Indiv1  1 1 1 0 -9
Indiv1  -9  -9  -9  -9
Indiv1  -9  -9  -9  -9
Indiv2  0 1 0 -9 0
Indiv2  0 1 0 -9 0
Indiv2  0 0 0 -9 0
Indiv2  1 0 0 -9 0
Indiv3  0 1 0 0 0
Indiv3  1 1 0 0 0
Indiv3  -9 1 1 -9 -9
Indiv3  -9 1 1 -9 -9
```

As with all scripts in this directory, options can be seen by calling the program with '-h':
```
./polyVCFtoStructure.py -h
Coming soon
```

