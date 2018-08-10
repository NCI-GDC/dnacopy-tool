library(data.table) 
library(dplyr)

setwd("~/SCRATCH/snp6/")

# required metadata files
cn.meta.hg19.file = "GenomeWideSNP_6.cn.na35.annot.bed"
cn.meta.hg38.file = "GenomeWideSNP_6.cn.na35.liftoverhg38.annot.bed"
snp.meta.hg19.file = "GenomeWideSNP_6.na35.annot.bed"
snp.meta.hg38.file = "GenomeWideSNP_6.na35.liftoverhg38.annot.bed"
nocnv.hg19.file = "nocnv.region.hg19.bed"
nocnv.hg38.file = "nocnv.region.hg38.bed"
par.hg38.file = "par.grch38.txt"

# read hg19 probe metadata 
colnames = c("chr", "start", "end", "strand", "ID") 
cn.meta.hg19 = fread(cn.meta.hg19.file, h=F, col.names = colnames, na.strings="---", )
snp.meta.hg19 = fread(snp.meta.hg19.file, h=F, col.names = colnames, na.strings="---", )
cn.meta.hg19$type = "CN"
snp.meta.hg19$type = "SNP"
meta.hg19 = rbind(cn.meta.hg19, snp.meta.hg19)
nocnv.hg19 = fread(nocnv.hg19.file, h=F, col.names=c("chr", "start", "end", "name"))

# read par region
par.hg38 = fread(par.hg38.file, col.names=c("name", "chr", "start", "end"))

# find nocnv probes in hg19 contect
setkey(nocnv.hg19, chr, start, end)
ov.hg19 = foverlaps(meta.hg19[!is.na(start)], nocnv.hg19, type="any", which=F, nomatch=0L)
meta.hg19$nocnv.hg19 = meta.hg19$ID %in% ov.hg19$ID
nocnv.hg38 = fread(nocnv.hg38.file, h=F, col.names=c("chr", "start", "end", "name"))


setwd("~/SCRATCH/snp6/meta201802/")
# load BWA fastmap remapped coordinates
snp = readRDS("snp6.snp.remap.rds")
names(snp)[1] = "ID"
snp$type = "SNP"
cn = readRDS("snp6.cn.remap.rds")
names(cn)[1] = "ID"
cn$type = "CN"
meta.hg38 = rbind(snp[, c("chr", "start", "end", "strand", "ID", "type", "reason"), with=F], cn[, c("chr", "start", "end", "strand", "ID", "type", "reason"), with=F])

# find nocnv probles in hg38 context
setkey(nocnv.hg38, chr, start, end)
ov.hg38 = foverlaps(meta.hg38[!is.na(start)], nocnv.hg38, type="any", which=F, nomatch=0L)
meta.hg38$freqcnv = meta.hg38$ID %in% intersect(ov.hg38$ID, ov.hg19$ID)


par.hg38$chr = paste0("chr", par.hg38$chr)
setkey(par.hg38, chr, start, end)
ov.par = foverlaps(meta.hg38[!is.na(start)], par.hg38, type="any", which=F, nomatch=0L)
meta.hg38$par = meta.hg38$ID %in% ov.par$ID

fwrite(meta.hg38, "snp6.na35.remap.hg38.full.txt", col.names=T, row.names=F, sep="\t", quote=F)


# add position and subset
meta.hg38$pos = as.integer(floor((meta.hg38$start + meta.hg38$end)/2))
names(meta.hg38)[which(names(meta.hg38)=="ID")] = "probeid"
meta.hg38 = meta.hg38[reason=="perfect_match"]
meta.hg38 = meta.hg38[, c("probeid", "chr", "pos", "strand", "type", "freqcnv", "par"), with=F]

# save meta
# save(meta, file="snp6.na35.liftoverhg38.rda") 
fwrite(meta.hg38, "snp6.na35.remap.hg38.txt", col.names=T, row.names=F, sep="\t", quote=F)

# 1880794 probes total
# in the previous liftover method, we have 1877330 that passed QC; and with remap method 1865120
# the intersect is 1863442. However, we recovered 1678 probes mainly due to probes map to a different 
# chromosome or A/B allele switch but still valid proble; we removed 13888 probes that mainly because 
# the same sequence now have perfect mapping to multiple locations, or both reference probe and alternative
# probe have perfect mappings to different locations. Together, that's close to 0.8% of the probes. 




