### MICROPEPTIDES ###
library(rtracklayer)
library(GenomicFeatures)
library(plyr)
library(data.table)

# get annotated ORFs shorter than 300nt

gtf_annot <- import("/Volumes/USELESS/DATA/genomes/GTF/Homo_sapiens.GRCh38.79.chr.gtf", format = "gtf")

txdb <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Homo_sapiens.GRCh38.79.chr.gtf", format = "gtf")

txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
## get transcript with longest cds per gene
#txlen <- arrange(txLengths, gene_id, desc(cds_len))
#txlen <- txlen[!duplicated(txlen$gene_id),]
#rownames(txlen) <- txlen$tx_name
#txlen <- txlen[order(rownames(txlen)),]

cds <- cdsBy(txdb, by="tx", use.names=TRUE)
exons <- exonsBy(txdb, by="tx", use.names=TRUE)

protein_coding <- unique(gtf_annot[gtf_annot$gene_biotype == "protein_coding"]$transcript_id)
#exons_pt <- exons[names(exons) %in% protein_coding]
#cds_pt <- cds[names(cds) %in% protein_coding]
short_cds <- txLengths[txLengths$cds_len <= 300 & txLengths$cds_len > 0,]$tx_name
sORFs <- cds[names(cds) %in% short_cds]
sORFs <- sORFs[names(sORFs) %in% protein_coding]

# remove overlapping annotations?
#txLengths[txLengths$tx_name %in% names(sORFs),]

# check for ribo-seq/RNA support in different libraries
ribo_FPKM <- get(load(file = "/Volumes/USELESS/DATA/MICROPEPTIDES/Ribo_FPKM.Rsave"))
rna_FPKM <- get(load(file = "/Volumes/USELESS/DATA/MICROPEPTIDES/RNA_FPKM.Rsave"))
rm(FPKM)

## TE - match ribo and rna samples



# sORFs with ribo-seq support
sORF_ribo <- ribo_FPKM[rownames(ribo_FPKM) %in% names(sORFs),]
sum(apply(sORF_ribo[,2:104], 1, function(x){any(x > 1)})) # 8197
sum(apply(sORF_ribo[,2:104], 1, function(x){all(x > 1)})) # 1236


###### lincRNA !!!!!!!!!!! remove all those overlapping pseudogenes etc. !!!!!!!!!!!!!!
lincRNA <- unique(gtf_annot[gtf_annot$gene_biotype == "lincRNA"]$transcript_id)
lincRNA_ribo <- ribo_FPKM[ribo_FPKM$tx %in% lincRNA,]

sum(apply(lincRNA_ribo[,2:104], 1, function(x){any(x > 1)})) # 1695
sum(apply(lincRNA_ribo[,2:104], 1, function(x){all(x > 1)})) # 6
# "ENST00000438002" "ENST00000453554" "ENST00000573479" "ENST00000602890" "ENST00000604849" "ENST00000617702"

lincRNA_rna <- rna_FPKM[rna_FPKM$tx %in% lincRNA,]


#################################################### TE
# ribo_TE <- ribo_FPKM
# ribo_TE <- ribo_TE[, -grep("^fritsch",colnames(ribo_TE)), with=FALSE]
# ribo_TE <- ribo_TE[, -grep("^ingolia",colnames(ribo_TE)), with=FALSE]
# ribo_TE <- ribo_TE[, -grep("^lee",colnames(ribo_TE)), with=FALSE]
# ribo_TE <- ribo_TE[, -grep("^liu",colnames(ribo_TE)), with=FALSE]
# ribo_TE <- ribo_TE[, -grep("^yoon",colnames(ribo_TE)), with=FALSE]
# rownames(ribo_TE) <- ribo_FPKM$tx
# ribo_TE$tx <- NULL
# 
# rna_TE <- rna_FPKM[, grep("^andreev",colnames(rna_FPKM)), with=FALSE]
# rna_TE <- cbind(rna_TE, rna_FPKM[, grep("^gonzalez",colnames(rna_FPKM)), with=FALSE])
# rna_TE <- cbind(rna_TE, rna_FPKM[,c("guo2010_SRR057513", "guo2010_SRR057514", "guo2010_SRR057518", "guo2010_SRR057519",
#                                     "guo2010_SRR057523", "guo2010_SRR057524", "guo2010_SRR057527", "guo2010_SRR057530",
#                                     "guo2010_SRR057533", "guo2010_SRR065776", "guo2010_SRR065777", "guo2010_SRR065781",
#                                     "guo2010_SRR065782")])
# rna_TE <- cbind(rna_TE, rna_FPKM[, grep("^hsieh",colnames(rna_FPKM)), with=FALSE])
# rna_TE <- cbind(rna_TE, rna_FPKM[,c("rutkowski2015_SRR1523653", "rutkowski2015_SRR1523653", "rutkowski2015_SRR1523653",
#                                     "rutkowski2015_SRR1523667", "rutkowski2015_SRR1523667", "rutkowski2015_SRR1523667")])
# rna_TE <- cbind(rna_TE, rna_FPKM[, grep("^sidrauski",colnames(rna_FPKM)), with=FALSE])
# rna_TE <- cbind(rna_TE, rna_FPKM[,c("stern2012_SRR592967", "stern2012_SRR592968", "stern2012_SRR592968",
#                                     "stern2012_SRR592966", "stern2012_SRR592968", "stern2012_SRR592966",
#                                     "stern2012_SRR592968", "stern2012_SRR592966", "stern2012_SRR592968",
#                                     "stern2012_SRR592966")])
# rna_TE <- cbind(rna_TE, rna_FPKM[, grep("^stumpf",colnames(rna_FPKM)), with=FALSE])
# rna_TE <- cbind(rna_TE, rna_FPKM[,c("subtelny2014_SRR1039860")])
# rownames(rna_TE) <- rna_FPKM$tx
# 
# 
# TE <- ribo_TE / rna_TE
# #TE_log2 <- log2(ribo_TE / rna_TE)
# save(TE, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/TE_FPKM.Rsave")

load(file = "/Volumes/USELESS/DATA/MICROPEPTIDES/TE_FPKM.Rsave")

####################################################


subtelny <- readGAlignments("/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/human/Subtelny__2014.Human.HEK293T.RPF.GRCh38.SRR1039861.bam")
#subtelny_exons <- subsetByOverlaps(subtelny, exons)
fo <- findOverlaps(subtelny, exons)
dark_matter <- subtelny[-from(fo)]

sort(table(start(dark_matter[seqnames(dark_matter) == "chr3"])), decreasing = TRUE)

rRNA_annot <- import("/Volumes/USELESS/DATA/genomes/GTF/human_rRNA.gtf", format = "gtf")




########## sORFs which do not have overlaps with other ORFs
ensembl_short_noncoding <- ensembl_short_noncoding[-unique(to(findOverlaps(ensembl_coding, ensembl_short_noncoding)))]

sORFs <- exons[names(exons) %in% short_cds]
sORFs <- sORFs[names(sORFs) %in% protein_coding]
#not_sORFs <- exons[names(exons)%in% protein_coding[!(protein_coding %in% names(sORFs))]]
not_sORFs <- exons[names(exons)%in% names(exons)[!(names(exons) %in% names(sORFs))]]

# find overlaps in sORFs and not_sORFs

unique_sORFs <- sORFs[-unique(from(findOverlaps(sORFs, not_sORFs)))] # multiple isoforms

# find sORFs without isoforms
txlen_sORFs <- txLengths[txLengths$tx_name %in% names(unique_sORFs),]
txlen_sORFs <- arrange(txlen_sORFs, gene_id, desc(cds_len))

ug <- table(txlen_sORFs$gene_id)
ug <- names(ug[ug == 1])
ut <- txlen_sORFs[txlen_sORFs$gene_id %in% ug,]$tx_name

unique_sORFs_noiso <- unique_sORFs[names(unique_sORFs) %in% ut]
#as.character(runValue(seqnames(unique_sORFs_noiso))) #chromosomes
sORFs_FPKM <- ribo_FPKM[ribo_FPKM$tx %in% ut]

################################################# on CCDS ################################

ccds <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/CCDS.bed")
names(ccds) <- ccds$name
ccds_len <- sapply(ccds$blocks, function(x){sum(width(x))}) # sum(ccds_len < 300) = 858
names(ccds_len) <- names(ccds)
short_ccds <- names(ccds_len[ccds_len < 300])
short_ccds <- ccds[names(ccds) %in% short_ccds]
n <- names(short_ccds)
short_ccds <- blocks(short_ccds)
names(short_ccds) <- n
short_ccds <- short_ccds[order(names(short_ccds))]

short_cds <- cds[unique(from(findOverlaps(cds, short_ccds)))]
txlen_short_cds <- txLengths[txLengths$tx_name %in% names(short_cds),]
# length(unique(txlen_short_cds$gene_id)) # 768

########### calculate FPKM on short_CCDS coordinates
# load each wiggle (bw) track, subset short_CCDS

wig_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/tracks"
wigs <- list.files(path = wig_path)

### import (separated by length?)
# andreev2015           "Andreev.+"
# fritsch2012           "Fritsch.+"
# gonzalez2014          "Gonzalez.+"
# guo2010               "Guo.+"
# hsieh2012             "Hsieh.+"
# ingolia2012           "Ingolia_NT_2012.+"
# ingolia2014           "Ingolia_NT_2014.+"
# lee2012               "Lee.+"
# liu2013               "Liu.+"
# rutkowski2015         "Rutkowski.+"
# sidrauski2015         "Sidrauski.+"
# stern2012             "Stern.+"
# stumpf2013            "Stumpf.+"
# subtelny2014          "Subtelny.+"
# yoon2014              "Yoon.+"

patterns <- c("Andreev.+", "Fritsch.+", "Gonzalez.+", "Guo.+", "Hsieh.+", "Ingolia_NT_2012.+", "Ingolia_NT_2014.+",
              "Lee.+", "Liu.+", "Rutkowski.+", "Sidrauski.+", "Stern.+", "Stumpf.+", "Subtelny.+", "Yoon.+")

prefixes <- c("andreev2015_", "fritsch2012_", "gonzalez2014_", "guo2010_", "hsieh2012_", "ingolia2012_",
              "ingolia2014_", "lee2012_", "liu2013_", "rutkowski2015_", "sidrauski2015_", "stern2012_",
              "stumpf2013_", "subtelny2014_", "yoon2014_")

short_ccds_len <- ccds_len[names(ccds_len) %in% names(short_ccds)]
short_ccds_len <- short_ccds_len[order(names(short_ccds_len))]


short_ccds_fpkm <- data.table(ccds_name = names(short_ccds), ccds_len = short_ccds_len)



for (i in 1:length(patterns)) {
  pattern <- patterns[i]
  prefix <- prefixes[i]
  
  libs <- list.files(path = wig_path, pattern = pattern)
  
  srr <- sapply(libs, function(x){strsplit(x, "[.]", fixed=F)})
  srr <- sapply(srr, function(x){x[length(x)-1]})
  srr <- unique(sapply(srr, function(x){substr(x, 1, nchar(x)-8)}))
  
  for (j in 1:length(srr)) {
    run <- srr[j]
    
    # fw
    pattern_fw <- paste0(pattern, run, "-forward.wig")
    wig_fw <- list.files(path = wig_path, pattern = pattern_fw)
    wig_fw <- import.wig(c(file.path(wig_path, wig_fw)))
    strand(wig_fw) <- rep("+", length(wig_fw))
    # rv
    pattern_rv <- paste0(pattern, run, "-reverse.wig")
    wig_rv <- list.files(path = wig_path, pattern = pattern_rv)
    wig_rv <- import.wig(c(file.path(wig_path, wig_rv)))
    strand(wig_rv) <- rep("-", length(wig_rv))
    # join
    wig <- c(wig_fw, wig_rv)
    seqlevels(wig)[which(seqlevels(wig) == "X")] <- "chrX"
    seqlevels(wig)[which(seqlevels(wig) == "Y")] <- "chrY"
    
    total_reads <- sum(wig$score)
    
    # map to short_CCDS
    wig_tx <- mapToTranscripts(wig, short_ccds, ignore.strand = FALSE)
    mcols(wig_tx) <- wig[wig_tx$xHits]$score
    
    # change to data.table, sum per CCDS
    wig_tx <- as.data.frame(wig_tx)
    setDT(wig_tx)
    
    # FPKM
    read_sum <- wig_tx[, .(sum = sum(X)), by = "seqnames"]
    colnames(read_sum) <- c("ccds_name", paste0(prefix, run))
    #cds_RNA_FPKM <- (cds_RNA / cds_len) * (10^9 / total_reads)
    short_ccds_fpkm <- merge(short_ccds_fpkm, read_sum, by = "ccds_name", all.x = TRUE)
    short_ccds_fpkm[[paste0(prefix,run)]] <- (short_ccds_fpkm[[paste0(prefix,run)]] / short_ccds_fpkm$ccds_len) * (10^9 / total_reads)
    
  }
  
}


short_ccds_fpkm[is.na(short_ccds_fpkm)] <- 0

save(short_ccds_fpkm, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/short_CCDS_FPKM.Rsave")

######################################################
############ openprot FPKM ###########################
openprot <- import(paste0("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/OpenProt/",
                          "human-openprot-r1_0-refprots+altprots+isoforms_min_1_pep-grch38.83+grch38.p7.bed/",
                          "human-openprot-r1_0-refprots+altprots+isoforms_min_1_pep-grch38.83+grch38.p7.bed"), format="bed")
n <- openprot$name
openprot <- blocks(openprot)
names(openprot) <- n
## short ones only
short_openprot <- openprot[sum(width(openprot)) < 300]

short_openprot_fpkm <- data.table(openprot_name = names(short_openprot), openprot_len = sum(width(short_openprot)))

for (i in 1:length(patterns)) {
  pattern <- patterns[i]
  prefix <- prefixes[i]
  
  libs <- list.files(path = wig_path, pattern = pattern)
  
  srr <- sapply(libs, function(x){strsplit(x, "[.]", fixed=F)})
  srr <- sapply(srr, function(x){x[length(x)-1]})
  srr <- unique(sapply(srr, function(x){substr(x, 1, nchar(x)-8)}))
  
  for (j in 1:length(srr)) {
    run <- srr[j]
    
    # fw
    pattern_fw <- paste0(pattern, run, "-forward.wig")
    wig_fw <- list.files(path = wig_path, pattern = pattern_fw)
    wig_fw <- import.wig(c(file.path(wig_path, wig_fw)))
    strand(wig_fw) <- rep("+", length(wig_fw))
    # rv
    pattern_rv <- paste0(pattern, run, "-reverse.wig")
    wig_rv <- list.files(path = wig_path, pattern = pattern_rv)
    wig_rv <- import.wig(c(file.path(wig_path, wig_rv)))
    strand(wig_rv) <- rep("-", length(wig_rv))
    # join
    wig <- c(wig_fw, wig_rv)
    seqlevels(wig)[which(seqlevels(wig) == "X")] <- "chrX"
    seqlevels(wig)[which(seqlevels(wig) == "Y")] <- "chrY"
    
    total_reads <- sum(wig$score)
    
    # map to short_CCDS
    wig_tx <- mapToTranscripts(wig, short_openprot, ignore.strand = FALSE)
    mcols(wig_tx) <- wig[wig_tx$xHits]$score
    
    # change to data.table, sum per CCDS
    wig_tx <- as.data.frame(wig_tx)
    setDT(wig_tx)
    
    # FPKM
    read_sum <- wig_tx[, .(sum = sum(X)), by = "seqnames"]
    colnames(read_sum) <- c("openprot_name", paste0(prefix, run))
    #cds_RNA_FPKM <- (cds_RNA / cds_len) * (10^9 / total_reads)
    short_openprot_fpkm <- merge(short_openprot_fpkm, read_sum, by = "openprot_name", all.x = TRUE)
    short_openprot_fpkm[[paste0(prefix,run)]] <- (short_openprot_fpkm[[paste0(prefix,run)]] / short_openprot_fpkm$openprot_len) * (10^9 / total_reads)
    
    print(run)
    
  }
  
}


short_openprot_fpkm[is.na(short_openprot_fpkm)] <- 0

save(short_openprot_fpkm, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/short_OpenProt_FPKM.Rsave")






