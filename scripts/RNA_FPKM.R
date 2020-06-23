## RNA-seq FPKM (for fantom_human)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(plyr)

### get transcript database
txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/Homo_sapiens.GRCh38.79.chr.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
exons <- exons[order(names(exons))]

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

## BAM files
libs_path <- "/export/valenfs/projects/fantom6/fancode/data/BAM/RNA-seq"
libs <- list.files(path = libs_path)
libs <- libs[-c(45,46)]

FPKM <- data.table(tx = sort(names(exons)))

for (i in 1:length(libs)) {
  ## load file
  RNAseq <- readGAlignments(libs[i])
  read_number_mRNA <- length(RNAseq)
  exons_RNA <- countOverlaps(exons, RNAseq)
  exons_RNA <- exons_RNA[order(names(exons_RNA))]
  exons_len <- txLengths[rownames(txLengths) %in% names(exons_RNA),]
  exons_len <- exons_len[order(rownames(exons_len)),]$tx_len
  RNA_FPKM <- (exons_RNA / exons_len) * (10^9 / read_number_mRNA)
  
  n <- unlist(strsplit(libs[i], split = "[.]"))
  lib <- paste0(n[1], "_", n[length(n)-1])
  FPKM[[lib]] <- RNA_FPKM
}

save(FPKM, file = "~/RNA_FPKM.Rsave")
save(FPKM, file = "/export/valenfs/projects/fantom6/fancode/data/BAM/RNA_FPKM.Rsave")

