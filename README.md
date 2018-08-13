# snp6cbs: Segment TCGA SNP6 tangent normalized data
- Version:	2.1
- Authors:	Zhenyu Zhang (zhenyuz.cdis@gmail.com)
- Description:	The package converts TCGA SNP6 tangent normalized data (produced previously by TCGA using BirdSuite package) into copy number segmentation files
- License: Apache 2.0

## Dependencies:	
- R (>= 3.2.0)
- tools
- parallel
- data.table (>= 1.9.6)
- DNAcopy (>= 1.44.0)
- futile.logger (>= 1.4.3)
- optparser (>= 1.4.4)


## Metadata Requirement
Probeset_Description File: a six column TSV file containing the following information
- probeid:	SNP6 probeset ID
- chr:	chromosome
- pos:	coordinate of the center of probe
- strand: strand
- type:	SNP6 or CN proble
- freqcnv: probe in frequent copy number variation regions in germline (boolean)
- par: probe in Pseudo-Autosomal Regions (PAR) (boolean)

Two example probe description files in GRCh38 are provided in ./meta folder. 
- snp6.na35.liftover.hg38.txt.gz is generated from coordinate liftover (default for command-line mode)
- snp6.na35.remap.hg38.subset.txt.gz is generated from probe sequence remap (recommended)

## Usage
Two modes of running are supported
In command-line running mode, the script SNP6CBS.r process one input file each time
```
Rscript SNP6CBS.r --file=tangent_normlized_file --out1=output_allcnv_file --out2=output_nocnv_file --meta=probeset_description_file --gender=gender --sample=sample_id --log=output_log_file 

```
In parallel mode, SNP6CBS.paralle.r will build a parallel cluster to process each sample. This mode requires a job description file, whose file name is hardcoded in the script as 
```
job_file = "../cbs.jobs.20180214.txt"
```
A sample job file is provided in ./example/cbs.jobs.20180214.txt. Please check file header for descriptions of each field.


