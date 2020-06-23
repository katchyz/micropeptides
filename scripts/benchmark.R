### benchmark - micropeptides
library(rtracklayer)
library(data.table)
library(ggplot2)
library(reshape2)
library(ORFik)

# load PRICE predictions, calculate fpkm on all predicted - see distribution
price_bed <- "/Users/kasia/Documents/PhD/scripts/micropeptides/PRICE/Price_1.0.2/price_bed/bed/price"
beds <- list.files(path = price_bed)
i <- 89 # ("stumpf2013_SRR970538.bed")
bed <- import(file.path(price_bed, beds[i]), format="bed")
price_orfs <- blocks(bed)

bed_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/bed_tracks_for_ORFik"
ribo <- list.files(path = bed_path)
i <- 95 # ("Stumpf_CR_2013.Human.HeLa.RPF.GRCh38.SRR970538.reads_merged.bed.gz")
rpf <- fread.bed(file.path(bed_path, ribo[i]))
seqlevels(rpf)[which(seqlevels(rpf) == "X")] <- "chrX"
seqlevels(rpf)[which(seqlevels(rpf) == "Y")] <- "chrY"
seqlevels(rpf)[which(seqlevels(rpf) == "MT")] <- "chrMT"

stumpf <- data.table(fpkm = fpkm(price_orfs, rpf), len = sum(width(price_orfs)))

# load sORFs.org predictions
sorfs_stumpf <- read.csv(file = "/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/sorfs_org/stumpf2013.txt", sep = "\t")

j <- 1
grl <- GRangesList()

for (i in 1:nrow(sorfs_stumpf)) {
  s <- sorfs_stumpf[i,]
  if (s$Spliced == "No") {
    chr <- paste0("chr", s$Chromosome)
    start <- s$Sorf.start
    end <- s$Sorf.end
    if (s$Strand == 1) {
      strand <- "+"
    } else if (s$Strand == -1) {
      strand <- "-"
    }
    
    gr <- GRanges(seqnames=Rle(chr, length(start)),
                  IRanges(start=start, end=end),
                  strand=strand)
    
  } else if (s$Spliced == "Yes") {
    chr <- paste0("chr", s$Chromosome)
    start <- as.numeric(unlist(strsplit(as.character(s$Spliced.start.parts), split = "_")))
    end <- as.numeric(unlist(strsplit(as.character(s$Spliced.stop.parts), split = "_")))
    if (s$Strand == 1) {
      strand <- "+"
    } else if (s$Strand == -1) {
      strand <- "-"
      temp <- start
      start <- rev(end)
      end <- rev(temp)
      rm(temp)
    }
    
    gr <- GRanges(seqnames=Rle(chr, length(start)),
            IRanges(start=start, end=end),
            strand=strand)
  }
  
  grl[[j]] <- gr
  j <- j+1
  
}

export.bed(grl, "/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/sorfs_org/stumpf2013.bed")


# short CCDS coordinates
ccds <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/CCDS.bed")
n <- ccds$name
ccds <- blocks(ccds)
names(ccds) <- n

ccds_len <- sum(width(ccds)) # sum(ccds_len < 300) = 858

short_ccds <- ccds[sum(width(ccds)) < 300]
short_ccds <- short_ccds[order(names(short_ccds))]
long_ccds <- ccds[sum(width(ccds)) >= 300]
long_ccds <- long_ccds[order(names(long_ccds))]

short_long <- findOverlaps(short_ccds, long_ccds)
short_ccds <- short_ccds[-unique(from(short_long))]
long_ccds <- long_ccds[-unique(to(short_long))]

# short CCDS FPKM
load(file = "/Volumes/USELESS/DATA/MICROPEPTIDES/short_CCDS_FPKM.Rsave") # short_ccds_fpkm
short_ccds_fpkm <- short_ccds_fpkm[short_ccds_fpkm$ccds_name %in% names(short_ccds),]

ccds_price_overlap <- data.table(ccds_name = names(short_ccds), ccds_len = sum(width(short_ccds)))
############ PRICE #############

price_bed <- "/Users/kasia/Documents/PhD/scripts/micropeptides/PRICE/Price_1.0.2/price_bed/bed/price"
beds <- list.files(path = price_bed)

for (i in 1:length(beds)) {
  lib <- substr(beds[i], 1, nchar(beds[i])-4)
  bed <- import(file.path(price_bed, beds[i]), format="bed")
  price_orfs <- blocks(bed)
  
  ccds_price <- findOverlaps(short_ccds, price_orfs)
  o <- rep("no", nrow(ccds_price_overlap))
  o[unique(from(ccds_price))] <- "yes"
  
  ccds_price_overlap[[paste0(lib, "_fpkm")]] <- short_ccds_fpkm[[lib]]
  ccds_price_overlap[[paste0(lib, "_overlap")]] <- o
  
  e <- mapply(function(x,y){all(c(start(x)==start(y), end(x)==end(y))) == TRUE}, short_ccds[from(ccds_price)], price_orfs[to(ccds_price)])
  
  ex <- rep("no_overlap", nrow(ccds_price_overlap))
  
  ex[unique(from(ccds_price)[which(e == F)])] <- "not_exact"
  ex[from(ccds_price)[which(e == T)]] <- "exact"
  
  ccds_price_overlap[[paste0(lib, "_exact")]] <- ex
  
}

save(ccds_price_overlap, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/CCDS_PRICE_OVERLAP.Rsave")





################ one by one lib ###################

for (i in 1:length(beds)) {
  lib <- substr(beds[i], 1, nchar(beds[i])-4)
  bed <- import(file.path(price_bed, beds[i]), format="bed")
  price_orfs <- blocks(bed)
  
  ccds_price_overlap <- data.table(ccds_name = names(short_ccds), ccds_len = short_ccds_len)
  
  ccds_price <- findOverlaps(short_ccds, price_orfs)
  o <- rep("no", nrow(ccds_price_overlap))
  o[unique(from(ccds_price))] <- "yes"
  
  ccds_price_overlap$fpkm <- short_ccds_fpkm[[lib]]
  ccds_price_overlap$overlap <- o
  
  e <- mapply(function(x,y){all(c(start(x)==start(y), end(x)==end(y))) == TRUE}, short_ccds[from(ccds_price)], price_orfs[to(ccds_price)])
  
  ex <- rep("not_found", nrow(ccds_price_overlap))
  
  ex[from(ccds_price)[which(e == T)]] <- "exact"
  ex[unique(from(ccds_price)[which(e == F)])] <- "not_exact"
  
  ccds_price_overlap$exact <- ex
  
  
  cols <- c("exact" = "red", "not_exact" = "#E69F00", "not_found" = "#999999")
  ggplot(ccds_price_overlap, aes(x=ccds_len, y=fpkm, colour=exact)) + geom_point() + scale_y_log10() + ggtitle(lib) +
    scale_color_manual(values=cols)
  
  ggsave(file=paste0("/Volumes/USELESS/DATA/MICROPEPTIDES/benchmark/PRICE/", lib, ".png"))
  
}


########## assess specificity (out-of-frame ORFs?) ########

dupa <- findOverlaps(short_ccds, short_ccds)
dupa <-dupa[!from(dupa) == to(dupa)]
unique_short_ccds <- short_ccds[-sort(unique(from(dupa)))]


################ one by one lib - unique CCDS ###################

for (i in 1:length(beds)) {
  lib <- substr(beds[i], 1, nchar(beds[i])-4)
  bed <- import(file.path(price_bed, beds[i]), format="bed")
  price_orfs <- blocks(bed)
  
  ccds_price_overlap <- data.table(ccds_name = names(short_ccds), ccds_len = sum(width(short_ccds)))
  
  ccds_price <- findOverlaps(short_ccds, price_orfs)
  o <- rep("no", nrow(ccds_price_overlap))
  o[unique(from(ccds_price))] <- "yes"
  
  ccds_price_overlap$fpkm <- short_ccds_fpkm[short_ccds_fpkm$ccds_name %in% names(short_ccds),][[lib]]
  ccds_price_overlap$overlap <- o
  
  e <- mapply(function(x,y){all(c(start(x)==start(y), end(x)==end(y))) == TRUE}, short_ccds[from(ccds_price)], price_orfs[to(ccds_price)])
  
  ex <- rep("not_found", nrow(ccds_price_overlap))
  
  ex[from(ccds_price)[which(e == T)]] <- "exact"
  ex[unique(from(ccds_price)[which(e == F)])] <- "not_exact"
  
  ccds_price_overlap$exact <- ex
  
  cols <- c("exact" = "red", "not_exact" = "#E69F00", "not_found" = "#999999")
  ggplot(ccds_price_overlap, aes(x=ccds_len, y=fpkm, colour=exact)) + geom_point() + scale_y_log10() + ggtitle(lib) +
    scale_color_manual(values=cols)
  
  ggsave(file=paste0("/Volumes/USELESS/DATA/MICROPEPTIDES/benchmark/price_unique_ccds/", lib, ".png"))
  
}


##############################################################################################################
##############################################################################################################
########## assess specificity (out-of-frame ORFs?) ########

# make sure short CCDSs do not overlap with longer CCDSs
ccds <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/CCDS.bed")
n <- ccds$name
ccds <- blocks(ccds)
names(ccds) <- n

ccds_len <- sum(width(ccds)) # sum(ccds_len < 300) = 858

short_ccds <- ccds[sum(width(ccds)) < 300]
short_ccds <- short_ccds[order(names(short_ccds))]
long_ccds <- ccds[sum(width(ccds)) >= 300]
long_ccds <- long_ccds[order(names(long_ccds))]

short_long <- findOverlaps(short_ccds, long_ccds)
short_ccds <- short_ccds[-unique(from(short_long))]
long_ccds <- long_ccds[-unique(to(short_long))]


for (i in 1:length(beds)) {
  lib <- substr(beds[i], 1, nchar(beds[i])-4)
  bed <- import(file.path(price_bed, beds[i]), format="bed")
  price_orfs <- blocks(bed)
  
  ccds_price_overlap <- data.table(ccds_name = names(short_ccds), ccds_len = sum(width(short_ccds)))
  
  ccds_price <- findOverlaps(short_ccds, price_orfs)
  o <- rep("no", nrow(ccds_price_overlap))
  o[unique(from(ccds_price))] <- "yes"
  
  ccds_price_overlap$fpkm <- short_ccds_fpkm[short_ccds_fpkm$ccds_name %in% names(short_ccds),][[lib]]
  ccds_price_overlap$overlap <- o
  
  e <- mapply(function(x,y){all(c(start(x)==start(y), end(x)==end(y))) == TRUE}, short_ccds[from(ccds_price)], price_orfs[to(ccds_price)])
  
  ex <- rep("not_found", nrow(ccds_price_overlap))
  
  ex[unique(from(ccds_price)[which(e == F)])] <- "not_exact"
  ex[from(ccds_price)[which(e == T)]] <- "exact"
  
  ccds_price_overlap$exact <- ex
  
  ### get not_exact, check if they are in frame
  not_exact_ccds <- short_ccds[from(ccds_price)[e==F]]
  not_exact_price <- price_orfs[to(ccds_price)[e==F]]
  names(not_exact_price) <- names(not_exact_ccds)
  
  # 
  nec <- tile1(not_exact_ccds)
  nep <- tile1(not_exact_price)
  #nec$frame <- rep(c(1,2,3))
  #nep$frame <- rep(c(1,2,3))
  fo <- mapply(function(x,y){findOverlaps(x,y)}, nec, nep)
  
  # from(fo)[1] == 1 # extended 5'end
  # to(fo)[1] == 1 # truncated 5'end
  
  ccds_price_overlap$inframe <- ccds_price_overlap$exact
  ## from(fo)[1] - to(fo)[1] %% 3 == 0 # in-frame (else: out-of-frame)
  ## length(nec) > length (nep) # truncated
  off_frame <- names(fo)[sapply(fo, function(x){from(x[1]) - to(x[1])}) %% 3 != 0]
  
  trun <- names(nec)[sapply(nec, length) > sapply(nep, length)]
  ext <- names(nec)[sapply(nec, length) < sapply(nep, length)]
  ccds_price_overlap[ccds_price_overlap$ccds_name %in% trun,]$inframe <- "truncation"
  ccds_price_overlap[ccds_price_overlap$ccds_name %in% ext,]$inframe <- "extention"
  ccds_price_overlap[ccds_price_overlap$ccds_name %in% off_frame,]$inframe <- "out-of-frame"
  
  ccds_price_overlap[ccds_price_overlap$inframe == "not_exact"]$inframe <- "diff_exon"
  ccds_price_overlap[ccds_price_overlap$exact == "exact"]$inframe <- "exact"
  
  cols <- c("exact" = "red", "extention" = "orange", "truncation" = "yellow", "not_found" = "gray", "out-of-frame" = "black", "diff_exon" = "pink")
  ggplot(ccds_price_overlap, aes(x=ccds_len, y=fpkm, colour=inframe)) + geom_point() + scale_y_log10() +
    ggtitle(lib) + scale_color_manual(values=cols)
  
  ggsave(file=paste0("/Volumes/USELESS/DATA/MICROPEPTIDES/benchmark/price_frame_trunc_ext/", lib, ".png"))
  
}



#######################
CPO <- get(load(file = "/Volumes/USELESS/DATA/MICROPEPTIDES/CCDS_PRICE_OVERLAP.Rsave"))
CPO <- CPO[CPO$ccds_name %in% names(short_ccds),]

CPO_fpkm <- CPO[,grepl( "_fpkm$" , colnames(CPO) ), with=FALSE]
CPO_exact <- CPO[,grepl( "_exact$" , colnames(CPO) ), with=FALSE]
rownames(CPO_fpkm) <- CPO$ccds_name
rownames(CPO_exact) <- CPO$ccds_name

DT <- stack(CPO_fpkm)
DT$exact <- stack(CPO_exact)$values
colnames(DT) <- c("fpkm", "study", "exact")
DT$study <- sapply(as.character(DT$study), function(x){substr(x, 1, nchar(x)-5)})
DT[DT$exact == "no_overlap",]$exact <- "not_found"
DT$orf <- CPO$ccds_name

ggplot(DT, aes(x=study, y=fpkm, fill=exact)) + geom_boxplot() + scale_y_log10()

stern <- DT[grepl( "^stern" , DT$study ),]
stern$study <- sapply(as.character(stern$study), function(x){substr(x, 11, nchar(x))})
ggplot(stern, aes(x=study, y=fpkm, fill=exact)) + geom_boxplot() + scale_y_log10()

stern$treatment <- stern$study
stern[stern$study == "SRR592954" | stern$study == "SRR592955" | stern$study == "SRR592956" | stern$study == "SRR609197",]$treatment <- "CHX"
stern[stern$study == "SRR592957" | stern$study == "SRR592958",]$treatment <- "Harr"
stern[stern$study == "SRR592959" | stern$study == "SRR592960",]$treatment <- "LTM"
stern[stern$study == "SRR592961" | stern$study == "SRR592962",]$treatment <- "no_drug"
ggplot(stern, aes(x=treatment, y=fpkm, fill=exact)) + geom_boxplot() + scale_y_log10() + ggtitle("stern2012")

stern$drug <- stern$treatment
stern[stern$study == "SRR592954",]$drug <- "CHX_1"
stern[stern$study == "SRR592955",]$drug <- "CHX_2"
stern[stern$study == "SRR592956",]$drug <- "CHX_3"
stern[stern$study == "SRR609197",]$drug <- "CHX_4"
stern[stern$study == "SRR592957",]$drug <- "Harr_1"
stern[stern$study == "SRR592958",]$drug <- "Harr_2"
stern[stern$study == "SRR592959",]$drug <- "LTM_1"
stern[stern$study == "SRR592960",]$drug <- "LTM_2"
stern[stern$study == "SRR592961",]$drug <- "no_drug_1"
stern[stern$study == "SRR592962",]$drug <- "no_drug_2"
ggplot(stern, aes(x=drug, y=fpkm, fill=exact)) + geom_boxplot() + scale_y_log10() + ggtitle("stern2012")

setDT(stern)
stern_exact <- stern[, sum(exact == "exact") ,by = orf]
stern_exact$V2 <- stern[, sum(exact == "not_exact") ,by = orf]$V1
stern_exact$V3 <- stern[, sum(exact == "not_found") ,by = orf]$V1
colnames(stern_exact) <- c("orf", "exact", "not_exact", "not_found")

mean_fpkm_exact <- stern[stern$exact == "exact",][, mean(fpkm), by=orf]
colnames(mean_fpkm_exact) <- c("orf", "mean_fpkm_exact")
mean_fpkm_notexact <- stern[stern$exact == "not_exact",][, mean(fpkm), by=orf]
colnames(mean_fpkm_notexact) <- c("orf", "mean_fpkm_notexact")
mean_fpkm_notfound <- stern[stern$exact == "not_found",][, mean(fpkm), by=orf]
colnames(mean_fpkm_notfound) <- c("orf", "mean_fpkm_notfound")

stern_exact <- merge(stern_exact, mean_fpkm_exact, by = "orf", all.x = TRUE)
stern_exact <- merge(stern_exact, mean_fpkm_notexact, by = "orf", all.x = TRUE)
stern_exact <- merge(stern_exact, mean_fpkm_notfound, by = "orf", all.x = TRUE)

stern_exact$freq_exact <- stern_exact$exact / (stern_exact$exact + stern_exact$not_exact)
stern_exact$freq_found <- (10 - stern_exact$not_found) / 10
stern_exact[stern_exact$mean_fpkm_notfound == 0,]$freq_found <- NaN

stern_exact$mean_fpkm_found <- (stern_exact$mean_fpkm_exact * stern_exact$exact +
  stern_exact$mean_fpkm_notexact * stern_exact$not_exact) / (stern_exact$exact + stern_exact$not_exact)

ggplot(stern_exact, aes(x=freq_exact)) + geom_histogram(bins=30)
ggplot(stern_exact, aes(x=freq_found)) + geom_histogram(bins=30)

ggplot(stern_exact, aes(x=freq_exact, y=mean_fpkm_exact)) + geom_point() + scale_y_log10()
ggplot(stern_exact, aes(x=freq_found, y=mean_fpkm_found)) + geom_point() + scale_y_log10()


##################
## check if PRICE-predicted ORFs have in-frame stop codons

# price_orfs
price_bed <- "/Users/kasia/Documents/PhD/scripts/micropeptides/PRICE/Price_1.0.2/price_bed/bed/price"
beds <- list.files(path = price_bed)
i <- 94
bed <- import(file.path(price_bed, beds[i]), format="bed")
price_orfs <- blocks(bed)
seqlevels(price_orfs)[which(seqlevels(price_orfs) == "chrMT")] <- "chrM"

# fasta_dna
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
seqlengths(genome)
#genome$chr1 # same as genome[["chr1"]]
seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, price_orfs, as.character = T)
minus <- seq[which(unlist(runValue(strand(price_orfs)) == "-"))]
minus <- sapply(minus, function(x){rev(x)})
seq[which(unlist(runValue(strand(price_orfs)) == "-"))] <- minus
seq <- DNAStringSet(sapply(seq, function(x){paste0(x, sep="", collapse="")}))
seq_translated <- translate(seq)



# extract sequence, divide by 3, look for taa/tga/tag





