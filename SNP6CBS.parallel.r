#!/usr/bin/env Rscript

# SNP6CBS Parallel 
# version: 2.1
# author: Zhenyu Zhang
# date: 08/10/2018
# require: R > 3.20 
#	   		data.table
#	   		DNAcopy
#			futile.logger
#			tools
#			parallel

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

format.segment = function(segment, gender=NA, sampleid = "Sample") {
	# format segment output, modify column names, and remove chrY is gender is female	
	names(segment) = c("GDC_Aliquot", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
	if (gender == "female") {
		segment = segment[segment$Chromosome != "Y", ]
	}
	segment$GDC_Aliquot = sampleid
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

# runCBS read input file and output segment files and log files
# tangent_file: input tangent copy number file
# meta_file: input gzipped snp6 probeset metadata name
# allcnv_file: output copy number segmentation file, default = "allcnv.txt"
# nocnv_file: output germline masked copy number segmentation file, default = "nocnv.txt"
# gender: gender, default = "unknown"
# sampleid: sample ID appear in the first column of the output, default = "Sample"
# log_file: output log file, default = "log.txt"
runCBS = function(tangent_file, meta_file, allcnv_file = "allcnv.txt", nocnv_file = "nocnv.txt", gender = "unknown", log_file = "log.txt", sampleid = "Sample") {
	
	######################################
	# Check and Read Input Data
	######################################

	# if gender is not "male" or "female", it will be set as "unknown"
	gender = tolower(gender)
	if (! gender %in% c("male", "female")) gender = "unknown"

	# start to write log 
	flog.appender(appender.file(log_file))

	# stop if input tangent_file is not provided or do not exist
	if (is.null(tangent_file)) flog.error("no input tangent normalized file provided")
	if (! file.exists(tangent_file)) {
		flog.error("input tangent normalized file does not exist")
		stop("Please provide tangent normalized file\n", call.=FALSE)
	}
	# read TCGA level 2 tangent CNV file
	tangentcn = read.tangentcn(tangent_file)

	# stop if input tangent_file is not provided or do not exist
	if (is.null(meta_file)) flog.error("no probe description file provided")
	if (! file.exists(meta_file)) {
		flog.error("input probe description file does not exist")
		stop("Please provide probe description file\n", call.=FALSE)
	}
	# read SNP6 metadata
	meta = read.meta(meta_file)

	# log all input options
	flog.info("input tangent file: \t\t\t%s", tangent_file)
	flog.info("input probe description file: \t\t\t%s", meta_file)
	flog.info("gender: \t\t\t%s", gender)
	flog.info("regular cnv output: \t\t%s", allcnv_file)
	flog.info("germline masked cnv output: \t%s", nocnv_file)
	flog.info("sample id used: \t\t%s", sampleid)


	######################################
	# Run CBS
	######################################

	# merge level 2 data with probeset metadata
	data = annotate.tangentcn(tangentcn = tangentcn, meta = meta)

	# exclude PAR if gender is male or unknown
	# this is disabled after we switch to probe positions that were carefully remapped
	# if (gender %in% c("male", "unknown")) {
	# 	data = data[par == FALSE]
	# }

	# calculate, format and write segment
	flog.info("Calculating regular CNV segments")
	mode = "allcnv"
	segs = get.segment(data, mode = mode, sampleid = sampleid)
	formatted.segs = format.segment(segs, gender = gender, sampleid = sampleid)
	write.table(formatted.segs, allcnv_file, col.names=T, row.names=F, sep="\t", quote=F)
	flog.info("Regular CNV Segnment Completed")

	# calculate, format and write nocnv segment
	flog.info("Calculating masked CNV segments")
	mode = "nocnv"
	nocnv.segs = get.segment(data, mode = mode, sampleid = sampleid)
	formatted.nocnv.segs = format.segment(nocnv.segs, gender = gender, sampleid = sampleid)
	write.table(formatted.nocnv.segs, nocnv_file, col.names=T, row.names=F, sep="\t", quote=F)
	flog.info("Masked CNV Segnment Completed")
}


# a wrapper function to load runCBS for production. 
# runCBSWrapper include data download and cleanup 
runCBSWrapper = function(tangent_uuid, tangent_file, meta_file, allcnv_file = "allcnv.txt", nocnv_file = "nocnv.txt", gender = "unknown", log_file = "log.txt", sampleid = "Sample", token_file){
	# Validate UUID
	if(grepl("^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", tangent_uuid)){
		# Download data
		cmd = paste("gdc-client download -t", token_file, "-d ./ --no-related-files", tangent_uuid)
		system(cmd)
		# Run
		tangent_file = paste0("./", tangent_uuid, "/", tangent_file)
		runCBS(tangent_file = tangent_file, meta_file = meta_file, allcnv_file = allcnv_file, nocnv_file = nocnv_file, gender = gender, sampleid = sampleid, log_file = log_file)
		# Remove downloaded data
		cmd = paste("rm -rf", tangent_uuid)
		system(cmd)
	}
}

# a wrapper function to use paralle to run multiple runCBSWrapper jobs
runCBSParallel = function(job_file, num_parallel){
	jobs = fread(job_file)
	cl = makeCluster(num_parallel, type="FORK")
	run = parLapply(cl, 1:nrow(jobs), function(x) 
			runCBSWrapper(	tangent_uuid = jobs$tangent_uuid[x], 
							tangent_file = jobs$tangent_file[x], 
							meta_file = jobs$meta_file[x], 
							allcnv_file = jobs$allcnv_file[x], 
							nocnv_file = jobs$nocnv_file[x], 
							gender = jobs$gender[x], 
							log_file = jobs$log_file[x], 
							sampleid = jobs$sample[x], 
							token_file = jobs$token_file[x]))	
	stopCluster(cl)
}


# main program to run multiple jobs parallely
require(data.table)
require(DNAcopy)
require(futile.logger)
require(tools)
require(parallel)

options(scipen=999)
setwd("~/SCRATCH/dev/cbs3/production20180214")

num_cores = detectCores()
job_file = "../cbs.jobs.20180214.txt"
runCBSParallel(job_file = job_file, num_parallel = num_cores - 10)




