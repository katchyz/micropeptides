## protein db (micropeptides)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(ggplot2)
library(ORFik)

## high confidence set SmProt, human
# smprot <- read.csv(file = "/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/SmProt/SmProt_human_highConf.txt", header = T, sep = "\t")
# 
# smprot$start <- as.integer(levels(smprot$start))[smprot$start]
# smprot$end <- as.integer(levels(smprot$end))[smprot$end]
# smprot <- smprot[-which(is.na(smprot$start)),] # missing start/end
# temp <- smprot[which(smprot$start > smprot$end),] # start or end disagree
# temp_switch <- temp[which(((temp$start - temp$end) + 1) / 3 == temp$length),]
# smprot[smprot$X.SmProt_ID %in% temp_switch$X.SmProt_ID,]$start <- temp_switch$end
# smprot[smprot$X.SmProt_ID %in% temp_switch$X.SmProt_ID,]$end <- temp_switch$start
# temp_remove <- temp[which(((temp$start - temp$end) + 1) / 3 != temp$length),]
# smprot <- smprot[-which(smprot$X.SmProt_ID %in% temp_remove$X.SmProt_ID),]
# rm(temp, temp_switch, temp_remove)
# smprot <- smprot[smprot$strand == "+" | smprot$strand == "-",] # missing strand
# 
# smprot <- GRanges(smprot)
# seqlevels(smprot)[which(seqlevels(smprot) == "chrM")] <- "chrMT"

######################################################
## CCDS: check alternative isoforms/overlaps
# CCDS.clustered.bed: 18703 clusters (32544 total)
######################################################
ccds <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/CCDS/CCDS.clustered.bed", format="bed")
ccds <- unique(split(ccds, ccds$NA.)) # 18703


######################################################
## OpenProt, detected with at least one unique peptide

# prior to loading: sort and cluster by bedtools (strand-specific)
# sort -k1,1 -k2,2n test.bed > test.sorted.bed
# bedtools cluster -s -i test.sorted.bed > test.clustered.bed

# openprot_all.clustered.bed: 32907 clusters
# short_openprot.clustered.bed: 14510 clusters
#####################################################

openprot <- import(paste0("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/OpenProt/",
                   "human-openprot-r1_0-refprots+altprots+isoforms_min_1_pep-grch38.83+grch38.p7.bed/",
                   "human-openprot-r1_0-refprots+altprots+isoforms_min_1_pep-grch38.83+grch38.p7.bed"), format="bed")
n <- openprot$name
openprot <- blocks(openprot) # 229800
names(openprot) <- n
## short ones only
short_openprot <- openprot[sum(width(openprot)) <= 300]
short_openprot <- short_openprot[order(names(short_openprot))] # 42374
long_openprot <- openprot[sum(width(openprot)) > 300]
long_openprot <- long_openprot[order(names(long_openprot))] # 187426

# remove overlaps with long proteins
overlap_short_long <- findOverlaps(short_openprot, long_openprot)
short_openprot <- short_openprot[-unique(from(overlap_short_long))] # 16628
long_openprot <- long_openprot[-unique(to(overlap_short_long))] # 81640

# export.bed(short_openprot, "/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/OpenProt/short_openprot.bed")
# remove multiple annotations (ensembl / refseq)
openprot_all_clust <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/OpenProt/openprot_all.clustered.bed")
openprot_all_clust <- unique(split(openprot_all_clust, openprot_all_clust$NA.)) #32907 genes

# single isoform







overlap_short_short <- findOverlaps(short_openprot, short_openprot)
overlap_short_short <- overlap_short_short[from(overlap_short_short) != to(overlap_short_short)]

# cluster isoforms
#unique_short_openprot <- ORFik::uniqueGroups(short_openprot)
#over_uso <- findOverlaps(unique_short_openprot, unique_short_openprot)
#over_uso <- over_uso[from(over_uso) != to(over_uso)]

i <- 1
l <- list()
for (idx in unique(from(overlap_short_short))) {
  indices <- c(idx, to(overlap_short_short[from(overlap_short_short) == idx]))
  l[[i]] <- setNames(indices, names(short_openprot[indices]))
  i <- i + 1
}

l <- sapply(l, function(x){sort(x)})
l <- unique(l)

# subset ?
keep <- c()
for (i in 1:length(l)) {
  keep <- c(keep, sum(sapply(l, function(x){all(l[[i]] %in% x)})) == 1)
}
l <- l[keep]

# merge non-pairwise-overlapping isoforms
merge <- list()
for (i in 1:length(l)) {
  merge[[i]] <- which(sapply(l, function(x){any(l[[i]] %in% x)}))
}

to_keep <- l[sapply(merge, length) == 1] # single
to_merge <- l[sapply(merge, length) > 1]

merge2 <- list()
for (i in 1:length(to_merge)) {
  merge2[[i]] <- which(sapply(to_merge, function(x){any(to_merge[[i]] %in% x)}))
}
merge2 <- unique(merge2)
for (i in 1:length(merge2)) {
  to_keep[[length(to_keep)+1]] <- setNames(unique(sort(unlist(to_merge[merge2[[i]]]))), unique(names(sort(unlist(to_merge[merge2[[i]]])))))
}

l <- to_keep

# check if the annotations are equal - keep unique records with only necessery identifiers
isoforms <- sapply(l, function(x){ORFik::uniqueGroups(short_openprot[x])})

# change names of isoforms
names_isoforms <- make.names(paste0(rep("gene", length(isoforms)), as.character(1:length(isoforms))))
names(isoforms) <- names_isoforms

for (i in 1:length(isoforms)) {
  g <- names(isoforms[i])
  l <- length(isoforms[[i]])
  t <- make.names(paste0(rep("iso", l), as.character(1:l)))
  gt <- make.names(paste0(rep(g, l), rep("_", l), t))
  names(isoforms[[i]]) <- gt
}

grl <- isoforms[[1]]
for (j in 2:length(isoforms)) {
  grl <- c(grl, isoforms[[j]])
} #geneX_isoY

# save(isoforms, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/OpenProt_short_gene_iso.Rsave")
# save(grl, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/OpenProt_short_isoforms.Rsave")

# calculate fpkm
bed_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/bed_tracks_for_ORFik"
beds <- list.files(path = bed_path)

OpenProt_FPKM <- data.table(gene_iso_name = names(grl), gene_iso_len = sum(width(grl)))
OpenProt_FPKM$gene_name <- sapply(names(grl), function(x){unlist(strsplit(x, split="_"))[1]})

for (i in 1:length(beds)) {
  rpf <- fread.bed(file.path(bed_path, beds[i]))
  seqlevels(rpf)[which(seqlevels(rpf) == "X")] <- "chrX"
  seqlevels(rpf)[which(seqlevels(rpf) == "Y")] <- "chrY"
  seqlevels(rpf)[which(seqlevels(rpf) == "MT")] <- "chrMT"
  # name
  bed <- unlist(strsplit(beds[i], split = "[.]"))
  surname <- paste0(tolower(unlist(strsplit(bed[1], split = "_"))[1]), "_", unlist(strsplit(bed[1], split = "_"))[3])
  srr <- bed[length(bed)-3]
  # fpkm
  OpenProt_FPKM[[paste0(surname, "_", srr)]] <- fpkm(grl, rpf)
}

# save(OpenProt_FPKM, file="/Volumes/USELESS/DATA/MICROPEPTIDES/OpenProt-FPKM.Rsave")






iso <- list()
for (j in 1:length(l)) {
  gene_clust <- l[[j]]
  iso_clust <- setNames(rep(F, length(gene_clust)), names(gene_clust))
  iso_clust[1] <- T
  for (k in 2:length(gene_clust)) {
    iso_clust[k] <- all(short_openprot[[gene_clust[1]]] == short_openprot[[gene_clust[k]]])
  }
  iso[[j]] <- iso_clust
}



####### check with PRICE
price_bed <- "/Users/kasia/Documents/PhD/scripts/micropeptides/PRICE/Price_1.0.2/price_bed/bed/price"
beds <- list.files(path = price_bed)

prot_price_overlap <- data.table(prot_name = names(short_openprot), prot_len = sum(width(short_openprot)))
load(file = "/Volumes/USELESS/DATA/MICROPEPTIDES/short_OpenProt_FPKM.Rsave")

for (i in 1:length(beds)) {
  lib <- substr(beds[i], 1, nchar(beds[i])-4)
  bed <- import(file.path(price_bed, beds[i]), format="bed")
  price_orfs <- blocks(bed)
  
  prot_price <- findOverlaps(short_openprot, price_orfs)
  o <- rep("no", nrow(prot_price_overlap))
  o[unique(from(prot_price))] <- "yes"
  
  #prot_price_overlap[[paste0(lib, "_fpkm")]] <- short_openprot_fpkm[[lib]]
  prot_price_overlap$fpkm <- short_openprot_fpkm[[lib]]
  #prot_price_overlap[[paste0(lib, "_overlap")]] <- o
  
  e <- mapply(function(x,y){all(c(start(x)==start(y), end(x)==end(y))) == TRUE}, short_openprot[from(prot_price)], price_orfs[to(prot_price)])
  
  ex <- rep("not_found", nrow(prot_price_overlap))
  
  ex[unique(from(prot_price)[which(e == F)])] <- "not_exact"
  ex[from(prot_price)[which(e == T)]] <- "exact"
  
  #prot_price_overlap[[paste0(lib, "_exact")]] <- ex
  prot_price_overlap$exact <- ex
  
  cols <- c("exact" = "red", "not_exact" = "#E69F00", "not_found" = "#999999")
  ggplot(prot_price_overlap, aes(x=prot_len, y=fpkm, colour=exact)) + geom_point(size=0.1) + scale_y_log10() + ggtitle(lib) +
    scale_color_manual(values=cols)
  
  ggsave(file=paste0("/Volumes/USELESS/DATA/MICROPEPTIDES/benchmark/price_openprot/", lib, ".png"))
  
}
