#!/usr/bin/env Rscript

# SNP6CBS
# author: Zhenyu Zhang
# date: 09/01/2016
# require: R > 3.20 
#	   	data.table
#	   	DNAcopy
#		futile.logger
#		optparse
#		tools

read.meta = function(file) {
	# read SNP6 metadata information in rda/rdata/gz/txt/tsv/csv format
	# get file extension
	flog.info("Loading metadata: %s", file)

	# read and format
	options(warn = -1)
	# gunzip to open gz data. zcat has compatibility issue in OSX
	meta = fread(paste("gunzip -c", file))
	options(warn = 0)
	meta$strand = factor(meta$strand)
	meta$type = factor(meta$type)
	meta$chr = factor(meta$chr, levels = c(1:22, "X", "Y"))
	# return
	return(meta)
}

get.segment = function(data, mode = "allcnv", sampleid = "Sample", seed=12345678) {
        # load data that at least contains segmean, chr, pos, return copy number segments
        if (mode == "nocnv") {
        	# remove freqcnv probes and chromosome Y
        	data = data[freqcnv == FALSE & chr != "Y"]
        }
        CNA.object = CNA( genomdat = data$segmean, chrom = data$chr, maploc = as.numeric(data$pos), data.type = "logratio", sampleid = sampleid)
        CNA.smoothed = smooth.CNA(CNA.object)
        set.seed(seed)
        segs = segment(CNA.smoothed, nperm=10000, alpha=0.01)
        return(segs$output)
}

format.segment = function(segment, gender=NA) {
	# format segment output, modify column names, and remove chrY is gender is female	
	names(segment) = c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
	if (gender == "female") {
		segment = segment[segment$Chromosome != "Y", ]
	}
	return(segment)
}

read.tangentcn = function(file) {
	# read tangentcn data
	flog.info("Loading tangent normalized copy number input file: %s", file)
	options(warn = -1)
	data = fread(file, skip=1, col.names=c("probeid", "signal"))
	options(warn = 0)
	return(data)	
}

annotate.tangentcn = function(tangentcn, meta) {
	# merge tangentcn data with meta data, and add segmean
	flog.info("Annotating tangent data with probeset information from metadata")
	data = merge(meta, tangentcn, by="probeid")
	# remove probes with missing information
	data = data[(! is.na(data$strand)) & (! is.na(data$signal)), ]
	# calculate probeset segment mean 
	data$segmean = log2(data$signal) - 1
	return(data)
}

get.opt = function() {
	# define option_list
	option_list = list(
		make_option(c("-f", "--file"), type = "character", default = NULL, 
	              	help = "tangent copy number file name"),
		make_option(c("--out1"), type = "character", default = "allcnv.txt", 
	              	help = "output copy number segmentation file name [default = %default]"), 
		make_option(c("--out2"), type = "character", default = "nocnv.txt", 
	              	help = "output germline masked copy number segmentation file name [default = %default]"), 
	    make_option(c("-m", "--meta"), type = "character", default = "snp6.na35.liftoverhg38.txt.gz", 
	              	help = "gzipped snp6 probeset metadata name [default = %default]"), 
	    make_option(c("--gender"), type = "character", default = "unknown", 
	              	help = "gender of the patient [default = %default]\n\t\t\t\t      [allowed = male, female, unknown]"),
	    make_option(c("-s", "--sample"), type = "character", default = "Sample", 
	              	help = "Sample ID in output files [default = %default]"), 
	    make_option(c("-l", "--log"), type = "character", default = "stdout", 
	              	help = "log file [default = %default]")   
	)
	
	# parse 
	opt_parser = OptionParser(option_list = option_list);
	opt = parse_args(opt_parser);

	# if log output is not "stdout", write to log file
	if(opt$log != "stdout") {
		flog.appender(appender.file(opt$log))
	}

	# stop if input file is not provided
	if (is.null(opt$file)){
	  print_help(opt_parser)
	  flog.error("no input file provided")
	  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
	}

	# change gender to unknown if not male or female
	opt$gender = tolower(opt$gender)
	if (! opt$gender %in% c("male", "female")) {
		opt$gender = "unknown"
	}

	# display all input options
	flog.info("input: \t\t\t%s", opt$file)
	flog.info("gender: \t\t\t%s", opt$gender)
	flog.info("regular cnv output: \t\t%s", opt$out1)
	flog.info("germline masked cnv output: \t%s", opt$out2)
	flog.info("metadata used: \t\t%s", opt$meta)
	flog.info("sample id used: \t\t%s", opt$sample)

	# return
	return(opt)
}


# load libraries
require(data.table)
require(DNAcopy)
require(optparse)
require(futile.logger)
require(tools)

options(scipen=999)

# get options
opt = get.opt()

# if log output is not "stdout", write to log file
if(opt$log != "stdout") {
	flog.appender(appender.file(opt$log))
}

# check input existance
# this script assume input file is valid TCGA tangent copy number file and will not check file content
if (! file.exists(opt$file)) {
	flog.error("input file does not exist")
	stop("Please provide input file\n", call.=FALSE)
}

# check metadata existance and file md5
# this script assume the metadata has md5sum fce277e26e4d65d187e1ea9800628bb9
md5 = md5sum(opt$meta)
if (is.na(md5)) {
	flog.error("metadata does not exist")
	stop("Please provide metadata\n", call.=FALSE)
} else if (md5 != "fce277e26e4d65d187e1ea9800628bb9") {
	flog.warning("metadata md5sum does not match fce277e26e4d65d187e1ea9800628bb9")
}

# read SNP6 metadata
meta = read.meta(opt$meta)
# read TCGA level 2 tangent CNV file
tangentcn = read.tangentcn(opt$file)
# merge level 2 data with probeset metadata
data = annotate.tangentcn(tangentcn = tangentcn, meta = meta)

# exclude PAR if gender is male or unknown
if (opt$gender %in% c("male", "unknown")) {
	data = data[par == FALSE]
}

# calculate, format and write segment
flog.info("Calculating regular CNV segments")
mode = "allcnv"
segs = get.segment(data, mode = mode, sampleid = opt$sample)
formatted.segs = format.segment(segs, gender = opt$gender)
write.table(formatted.segs, opt$out1, col.names=T, row.names=F, sep="\t", quote=F)
flog.info("Regular CNV Segnment Completed")

# calculate, format and write nocnv segment
flog.info("Calculating masked CNV segments")
mode = "nocnv"
nocnv.segs = get.segment(data, mode = mode, sampleid = opt$sample)
formatted.nocnv.segs = format.segment(nocnv.segs, gender = opt$gender)
write.table(formatted.nocnv.segs, opt$out2, col.names=T, row.names=F, sep="\t", quote=F)
flog.info("Masked CNV Segnment Completed")


