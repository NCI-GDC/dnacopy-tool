# script to prepare for SNP6 metadata for GRCh38 segmentation

# LiftOver SNP6 probeset coordinates from hg19 to GRCh38 using CrossMap
# Step 1 (Optional): create header row if desired later to concartenate with BED files
# > echo -e "chrom\tstart\tend\tstrand\tID"  > header.txt

# Step 2: prepare input for CrossMap
# > zcat GenomeWideSNP_6.cn.na35.annot.csv.zip | tail -n+20 | sed 's/"//g' | awk -F"," '{print "chr"$2"\t"$3"\t"$4"\t"$5"\t"$1}' > GenomeWideSNP_6.cn.na35.annot.bed
# > zcat GenomeWideSNP_6.na35.annot.csv.zip | tail -n+20 | sed 's/"//g' | awk -F"," '{print "chr"$3"\t"$4"\t"$4"\t"$5"\t"$1}' > GenomeWideSNP_6.na35.annot.bed
# > tail -n+2 nocnv.region.hg19.txt | awk '{print "chr"$2"\t"$3"\t"$4"\t"$1}' | sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g' > nocnv.region.hg19.bed

# Step 3: liftover coordinates to hg38
# export input_chain_file="$HOME/CrossMap-0.2.2/data/hg19ToHg38.over.chain.gz"
# python $HOME/CrossMap/usr/local/bin/CrossMap.py bed $input_chain_file GenomeWideSNP_6.cn.na35.annot.bed GenomeWideSNP_6.cn.na35.liftoverhg38.annot.bed
# python $HOME/CrossMap/usr/local/bin/CrossMap.py bed $input_chain_file GenomeWideSNP_6.na35.annot.bed GenomeWideSNP_6.na35.liftoverhg38.annot.bed

# Step 4: Prepare SNP6 metadata using
#	 			Affymetrics metadata, 
#				hg38 liftOver metadata,
#         		Intervals of Pseudo-Autosomal Regions on hg38
#				Frequent variant regions of normal copy numbers from Broad
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
setwd("./metadata/")

# required metadata files
cn.meta.hg19.file = "GenomeWideSNP_6.cn.na35.annot.bed"
cn.meta.hg38.file = "GenomeWideSNP_6.cn.na35.liftoverhg38.annot.bed"
snp.meta.hg19.file = "GenomeWideSNP_6.na35.annot.bed"
snp.meta.hg38.file = "GenomeWideSNP_6.na35.liftoverhg38.annot.bed"
nocnv.hg19.file = "nocnv.region.hg19.bed"
nocnv.hg38.file = "nocnv.region.hg38.bed"
par.hg38.file = "par.grch38.txt"

# read hg19 probe metadata and CrossMap liftover coordinates
colnames = c("chr", "start", "end", "strand", "ID") 
cn.meta.hg19 = fread(cn.meta.hg19.file, h=F, col.names = colnames, na.strings="---", )
cn.meta.hg38 = fread(cn.meta.hg38.file, h=F, col.names = colnames, na.strings="---", )
snp.meta.hg19 = fread(snp.meta.hg19.file, h=F, col.names = colnames, na.strings="---", )
snp.meta.hg38 = fread(snp.meta.hg38.file, h=F, col.names = colnames, na.strings="---", )

# merge hg19/hg38, and Affy/LiftOver meta data together: 1878043 total
cn.meta = merge(cn.meta.hg38, cn.meta.hg19, by="ID", suffixes = c(".hg38", ".hg19"))
cn.meta$type = "CN"
snp.meta = merge(snp.meta.hg38, snp.meta.hg19, by="ID", suffixes = c(".hg38", ".hg19"))
snp.meta$type = "SNP"
meta = rbind(cn.meta, snp.meta)

# remove duplicated probes (cn probes mapped to multiple locations)
dupids = meta$ID[duplicated(meta$ID)]
meta = meta[! meta$ID %in% dupids, ]

# remove liftover coordinates not on major chromosomes:
chrs = c(1:22, "X", "Y")
meta = meta[chr.hg38 %in% paste0("chr", chrs), ]

# remove obvious probeset errors in hg38 by comparing the length of mapped probe length
# (a better but more complicated method is to map/liftover probeset sequences to hg38)
meta = meta[(end.hg38 - start.hg38) == (end.hg19 - start.hg19), ]

# remove liftover chromosomal changes (optional)
meta = meta[chr.hg19 == chr.hg38, ]

# make genomic ranges object for interval comparison
meta$chr = factor(gsub("chr", "", meta$chr.hg38), levels = chrs)
meta$strand.hg38[is.na(meta$strand.hg38)] = "*"
probe.hg38 = with(meta, GRanges(	seqname = chr, ranges = IRanges(start = start.hg38, end = end.hg38), 
							strand = strand.hg38, name = ID))
probe.hg19 = with(meta, GRanges(	seqname = chr, ranges = IRanges(start = start.hg19, end = end.hg19), 
							strand = "*", name = ID))

# read hg19 nocnv intervals, fix chromosome names, and make genomicranges object
nocnv.hg19 = fread(nocnv.hg19.file, h=F, col.names=c("chr", "start", "end", "name"))
nocnv.hg19$chr = gsub("chr", "", nocnv.hg19$chr)
nocnv.hg19 = with(nocnv.hg19, GRanges(	seqname = chr, ranges = IRanges(start = start, end = end), 
								strand = "*", name = name))

# annotate probeset by hg19 nocnv intervals
Hits = findOverlaps(probe.hg19 , nocnv.hg19, ignore.strand=T)
meta$freqcnv.hg19 = F
meta$freqcnv.hg19[queryHits(Hits)] = T


# read hg38 nocnv intervals, fix chromosome names, and make genomicranges object
nocnv.hg38 = fread(nocnv.hg38.file, h=F, col.names=c("chr", "start", "end", "name"))
nocnv.hg38$chr = gsub("chr", "", nocnv.hg38$chr)
nocnv.hg38 = with(nocnv.hg38, GRanges(	seqname = chr, ranges = IRanges(start = start, end = end), 
								strand = "*", name = name))

# annotate probeset by hg38 nocnv intervals
Hits = findOverlaps(probe.hg38 , nocnv.hg38, ignore.strand=T)
meta$freqcnv.hg38 = F
meta$freqcnv.hg38[queryHits(Hits)] = T

# merge nocnv probes
meta$freqcnv = meta$freqcnv.hg19 | meta$freqcnv.hg38

# read pseudo-autosomal region (PAR) intervals, fix chromosome names, and make genomicranges object
par.hg38 = fread(par.hg38.file, col.names=c("name", "chr", "start", "end"))
par.hg38 = with(par.hg38, GRanges(	seqname = chr, ranges = IRanges(start = start, end = end), 
								strand = "*", name = name))

# annotate probeset by nocnv intervals
Hits = findOverlaps(probe.hg38, par.hg38, ignore.strand=T)
meta$par = F
meta$par[queryHits(Hits)] = T

# add position
meta$pos = as.integer(floor((meta$start.hg38 + meta$end.hg38)/2))

# save metadata
# save(meta, file="snp6.na35.liftoverhg38.full.rda") 
write.table(meta, "snp6.na35.liftoverhg38.full.txt", col.names=T, row.names=F, sep="\t", quote=F)

# subsetting metadata
colnames = c("ID", "chr", "pos", "strand.hg38", "type", "freqcnv", "par")
meta = meta[, colnames, with=F]
colnames = c("probeid", "chr", "pos", "strand", "type", "freqcnv", "par")
names(meta) = colnames
meta$strand = factor(meta$strand)
meta$type = factor(meta$type)

# save meta
# save(meta, file="snp6.na35.liftoverhg38.rda") 
write.table(meta, "snp6.na35.liftoverhg38.txt", col.names=T, row.names=F, sep="\t", quote=F)


