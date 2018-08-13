# script to prepare for SNP6 metadata for GRCh38 segmentation

# Step 1: prepare input for CrossMap
# > zcat GenomeWideSNP_6.cn.na36.annot.csv.zip | tail -n+20 | sed 's/"//g' | awk -F"," '{print "chr"$2"\t"$3"\t"$4"\t"$5"\t"$1}' > GenomeWideSNP_6.cn.na36.annot.bed
# > zcat GenomeWideSNP_6.na36.annot.csv.zip | tail -n+20 | sed 's/"//g' | awk -F"," '{print "chr"$3"\t"$4"\t"$4"\t"$5"\t"$1}' > GenomeWideSNP_6.na36.annot.bed
# > tail -n+2 nocnv.region.hg19.txt | awk '{print "chr"$2"\t"$3"\t"$4"\t"$1}' | sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g' > nocnv.region.hg19.bed


# Step 2: Prepare SNP6 metadata using
#	 			Affymetrics metadata, 
#				hg38  metadata,
#         		Intervals of Pseudo-Autosomal Regions on hg38
#				hg38 liftOver frequent variant regions


require(data.table) 
require(GenomicRanges)

# if package not found
#
#	install.packages("data.table")
#
# 	source("http://bioconductor.org/biocLite.R")
#	biocLite("GenomicRanges")

# init
setwd("/mnt/SCRATCH/snp6/na36")

# required metadata files
cn.meta.hg38.file = "GenomeWideSNP_6.cn.na36.annot.bed"
snp.meta.hg38.file = "GenomeWideSNP_6.na36.annot.bed"
nocnv.hg38.file = "../nocnv.region.hg38.bed"
par.hg38.file = "../par.grch38.txt"

# read hg19 probe metadata and CrossMap liftover coordinates
colnames = c("chr", "start", "end", "strand", "ID") 
cn.meta.hg38 = fread(cn.meta.hg38.file, h=F, col.names = colnames, na.strings="---", )
snp.meta.hg38 = fread(snp.meta.hg38.file, h=F, col.names = colnames, na.strings="---", )

# merge SNP and CN meta data together: 1877768 total
cn.meta.hg38$type = "CN"
snp.meta.hg38$type = "SNP"
meta = rbind(cn.meta.hg38, snp.meta.hg38)

# remove liftover coordinates not on major chromosomes:
chrs = c(1:22, "X", "Y")
meta = meta[chr %in% paste0("chr", chrs), ]

# make genomic ranges object for interval comparison
meta$chr = factor(gsub("chr", "", meta$chr), levels = chrs)
probe = with(meta, GRanges(	seqname = chr, ranges = IRanges(start = start, end = end), 
							strand = strand, name = ID))

# read hg38 nocnv intervals, fix chromosome names, and make genomicranges object
nocnv.hg38 = fread(nocnv.hg38.file, h=F, col.names=c("chr", "start", "end", "name"))
nocnv.hg38$chr = gsub("chr", "", nocnv.hg38$chr)
nocnv.hg38 = with(nocnv.hg38, GRanges(	seqname = chr, ranges = IRanges(start = start, end = end), 
								strand = "*", name = name))

# annotate probeset by hg38 nocnv intervals
Hits = findOverlaps(probe , nocnv.hg38, ignore.strand=T)
meta$freqcnv = F
meta$freqcnv[queryHits(Hits)] = T

# read pseudo-autosomal region (PAR) intervals, fix chromosome names, and make genomicranges object
par.hg38 = fread(par.hg38.file, col.names=c("name", "chr", "start", "end"))
par.hg38 = with(par.hg38, GRanges(	seqname = chr, ranges = IRanges(start = start, end = end), 
								strand = "*", name = name))

# annotate probeset by nocnv intervals
Hits = findOverlaps(probe, par.hg38, ignore.strand=T)
meta$par = F
meta$par[queryHits(Hits)] = T

# add position
meta$pos = as.integer(floor((meta$start + meta$end)/2))


write.table(meta, "snp6.na36.full.txt", col.names=T, row.names=F, sep="\t", quote=F)

# subsetting metadata
colnames = c("ID", "chr", "pos", "strand", "type", "freqcnv", "par")
meta = meta[, colnames, with=F]
colnames = c("probeid", "chr", "pos", "strand", "type", "freqcnv", "par")
names(meta) = colnames
meta$strand = factor(meta$strand)
meta$type = factor(meta$type)

# save meta
# save(meta, file="snp6.na36.rda") 
write.table(meta, "snp6.na36.txt", col.names=T, row.names=F, sep="\t", quote=F)


