### Ribo_FPKM ###
library(rtracklayer)
library(GenomicFeatures)
library(plyr)
library(data.table)

txdb <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/Homo_sapiens.GRCh38.79.chr.gtf", format = "gtf")

exons <- exonsBy(txdb, by="tx", use.names=TRUE)
exons <- exons[order(names(exons))]

txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]


# check for ribo-seq/RNA support in different libraries
libs_path <- "/export/valenfs/projects/fantom6/fancode/data/fantom_human_tracks/granges"
libs <- list.files(path = libs_path)

FPKM <- data.table(tx = sort(names(exons)))

for (i in 1:length(libs)) {
  ## load file
  load(file = c(file.path(libs_path, libs[i])))
  grdt <- as.data.frame(gr)
  setDT(grdt)
  rm(gr)
  
  ## sum riboseq per tx
  riboseq <- grdt[, .(length = (.N), sum = sum(riboseq)), by=seqnames]
  footprint_number_total <- sum(riboseq$sum)
  
  ribo_FPKM <- (riboseq$sum / riboseq$length) * (10^9 / footprint_number_total)
  
  lib <- substr(libs[i], 1, nchar(libs[i])-6)
  FPKM[[lib]] <- ribo_FPKM
}

save(FPKM, file = "~/Ribo_FPKM.Rsave")
save(FPKM, file = "/export/valenfs/projects/fantom6/fancode/data/BAM/Ribo_FPKM.Rsave")

