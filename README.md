# snp6cbs

## Tool Description
Package:	snp6cbs
Title:	Segment TCGA SNP6 tangent normalized data
Version:	2.1
Authors:	Zhenyu Zhang (zhenyuz.cdis@gmail.com)
Description:	The package converts TCGA SNP6 tangent normalized data (produced previously by TCGA using BirdSuite package) into copy number segmentation files
Depends:	R (>= 3.2.0)
			tools
			parallel
			data.table (>= 1.9.6)
			DNAcopy (>= 1.44.0)
			futile.logger (>= 1.4.3)
			optparser (>= 1.4.4)
License: Apache 2.0

## Input Requirement
Probe Description File: a six column TSV file containing the following information
		probeid:	SNP6 probeset ID
		chr:	chromosome
		pos:	coordinate of the center of probe
		strand: strand
		type:	SNP6 or CN proble
		freqcnv: probe in frequent copy number variation regions in germline (boolean)
		par: probe in Pseudo-Autosomal Regions (PAR) (boolean)

