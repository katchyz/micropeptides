### compute features on short and long CCDSs
library(rtracklayer)
library(ORFik)
library(GenomicFeatures)
library(data.table)
library(ggplot2)

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

######
#short_ccds <- makeORFNames(short_ccds)
#long_ccds <- makeORFNames(long_ccds)

bed_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/bed_tracks_for_ORFik"
beds <- list.files(path = bed_path)

gtf <- makeTxDbFromGFF("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/CCDS.gtf", format = "gtf")

for (i in 1:length(beds)) {
  rpf <- fread.bed(file.path(bed_path, beds[i]))
  seqlevels(rpf)[which(seqlevels(rpf) == "X")] <- "chrX"
  seqlevels(rpf)[which(seqlevels(rpf) == "Y")] <- "chrY"
  minL <- min(rpf$score)
  maxL <- max(rpf$score)
  
  short_ccds_features <- computeFeatures(short_ccds, rpf, RNA = NULL, Gtf = gtf, riboStart = minL, riboStop = maxL)
  
  fpkm(short_ccds, rpf)
  entropy(short_ccds, rpf)
  floss(short_ccds, rpf, short_ccds, start = minL, end = maxL)
  #fractionLength(short_ccds, sum(width(short_ccds)))
  #disengagementScore(grl, RFP, GtfOrTx)
  #ribosomeReleaseScore(grl, RFP, GtfOrThreeUtrs, RNA = NULL)
  ribosomeStallingScore(short_ccds, rpf)
  orfScore(short_ccds, rpf)
  #distToCds(ORFs, fiveUTRs, cds = NULL, extension = NULL)
  #kozakSequenceScore(grl, faFile, species = "human", include.N = FALSE)
  #insideOutsideORF(grl, RFP, GtfOrTx)
  
  
}

?fread.bed # can read gzipped
?computeFeatures

dupa <- c(short_ccds, long_ccds)
dupa <- dupa[order(names(dupa))]

dupa_fpkm <- fpkm(dupa, rpf)
dupa_entropy <- entropy(dupa, rpf)
dupa_floss <- floss(dupa, rpf, dupa, start = minL, end = maxL)
dupa_rss <- ribosomeStallingScore(dupa, rpf)
#dupa <- removeMetaCols(dupa)
dupa_orfscore <- orfScore(dupa, rpf)

dt <- data.table(name = names(dupa), fpkm = dupa_fpkm, entropy = dupa_entropy, floss = dupa_floss, rss = dupa_rss,
                 orfscore = dupa_orfscore$ORFScore)
dt$len <- sum(width(dupa))
dt$length <- rep("long", nrow(dt))
dt[dt$len < 300]$length <- "short"

ggplot(dt, aes(x=fpkm, y=entropy, colour=length)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(dt, aes(x=fpkm, y=floss, colour=length)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(dt, aes(x=fpkm, y=rss, colour=length)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(dt, aes(x=fpkm, y=orfscore, colour=length)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()

ggplot(dt, aes(x=len, y=entropy, colour=length)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(dt, aes(x=len, y=floss, colour=length)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(dt, aes(x=len, y=rss, colour=length)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(dt, aes(x=len, y=orfscore, colour=length)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()

ggplot(dt, aes(x=len, y=fpkm, colour=length)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()


####################################################################################################################
############ get transcripts with single isoform (whole tx, not just CDS)
############ subset those with ribo coverage - check if all exons have coverage
############ plot tx_profile
############ calculate features on CDSs
##############################################################################
############ get off-frame ORF on single isoform genes
############ calculate features on off-frame ORFs
##############################################################################
############ plot features together for CDSs and off-frame ORFs

OpenProtClust <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/OpenProt/short_openprot.clustered.bed")
OpenProtClust <- unique(split(OpenProtClust, OpenProtClust$NA.)) # 14510

openprot_one_iso <- OpenProtClust[sapply(OpenProtClust, function(x){length(x) == 1})]
openprot_one_iso <- unlist(openprot_one_iso)

cds_ensembl <- import.bed("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/Ensembl/cds_ensembl.bed")
seqlevels(cds_ensembl)[which(seqlevels(cds_ensembl) == "X")] <- "chrX"
seqlevels(cds_ensembl)[which(seqlevels(cds_ensembl) == "Y")] <- "chrY"
seqlevels(cds_ensembl)[which(seqlevels(cds_ensembl) == "MT")] <- "chrMT"

one_iso_ensembl_openprot <- unique(cds_ensembl[cds_ensembl %in% openprot_one_iso])$name

####################################################################################################################
################################# clustering Ensembl - one/two isoforms ############################################
####################################################################################################################

EnsemblClust <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/Ensembl/tx_ensembl.clustered.bed")
#EnsemblClust <- unique(split(EnsemblClust, EnsemblClust$NA.)) # 20485
one_iso_cluster_numbers <- names(table(EnsemblClust$NA.)[table(EnsemblClust$NA.) == 1])
ensembl_one_iso <- EnsemblClust[EnsemblClust$NA. %in% one_iso_cluster_numbers]
#ensembl_one_iso <- EnsemblClust[sapply(EnsemblClust, function(x){length(x) == 1})] # 4584
#ensembl_one_iso <- unlist(ensembl_one_iso)
seqlevels(ensembl_one_iso)[which(seqlevels(ensembl_one_iso) == "X")] <- "chrX"
seqlevels(ensembl_one_iso)[which(seqlevels(ensembl_one_iso) == "Y")] <- "chrY"
seqlevels(ensembl_one_iso)[which(seqlevels(ensembl_one_iso) == "MT")] <- "chrMT"
names(ensembl_one_iso) <- ensembl_one_iso$name
# get only chromosomes (not scaffolds)
ensembl_one_iso <- ensembl_one_iso[grepl("^chr", as.character(seqnames(ensembl_one_iso)))]
ensembl_one_iso <- ensembl_one_iso[order(names(ensembl_one_iso))]
n <- ensembl_one_iso$name
ensembl_one_iso <- blocks(ensembl_one_iso)
names(ensembl_one_iso) <- n

# stumpf2013_SRR970538
bed_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/bed_tracks_for_ORFik"
beds <- list.files(path = bed_path)
i <- 95
rpf <- fread.bed(file.path(bed_path, beds[i]))
seqlevels(rpf)[which(seqlevels(rpf) == "X")] <- "chrX"
seqlevels(rpf)[which(seqlevels(rpf) == "Y")] <- "chrY"
seqlevels(rpf)[which(seqlevels(rpf) == "MT")] <- "chrMT"
ensembl_fpkm <- fpkm(ensembl_one_iso, rpf)
ensembl_fpkm <- data.table(tx = names(ensembl_one_iso), fpkm = ensembl_fpkm)

one_iso_tx <- ensembl_fpkm[ensembl_fpkm$fpkm > 1,]$tx # 753
# load gtf and txdb -> make CDSs
txdb <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Homo_sapiens.GRCh38.79.chr.gtf", format = "gtf")
cds <- cdsBy(txdb, by="tx", use.names=TRUE)
exons <- exonsBy(txdb, by="tx", use.names=TRUE)
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)

tx_len <- txLengths[txLengths$tx_name %in% one_iso_tx,]
one_iso_tx <- sort(tx_len[tx_len$utr5_len > 0,]$tx_name) # need to have leaders # 697
one_iso_cds <- cds[names(cds) %in% one_iso_tx]
one_iso_cds <- one_iso_cds[order(names(one_iso_cds))] # 697
one_iso_exons <- exons[names(exons) %in% one_iso_tx]
one_iso_exons <- one_iso_exons[order(names(one_iso_exons))] # 697

seqlevels(one_iso_cds)[which(seqlevels(one_iso_cds) == "X")] <- "chrX"
seqlevels(one_iso_cds)[which(seqlevels(one_iso_cds) == "Y")] <- "chrY"
seqlevels(one_iso_cds)[which(seqlevels(one_iso_cds) == "MT")] <- "chrMT"
seqlevels(one_iso_exons)[which(seqlevels(one_iso_exons) == "X")] <- "chrX"
seqlevels(one_iso_exons)[which(seqlevels(one_iso_exons) == "Y")] <- "chrY"
seqlevels(one_iso_exons)[which(seqlevels(one_iso_exons) == "MT")] <- "chrMT"

#one_iso_features <- computeFeatures(one_iso, rpf, RNA = NULL, Gtf = txdb, riboStart = min(rpf$score), riboStop = max(rpf$score))

one_iso_fpkm <- fpkm(one_iso_cds, rpf)
one_iso_entropy <- entropy(one_iso_cds, rpf)
one_iso_floss <- floss(one_iso_cds, rpf, one_iso_cds, start = min(rpf$score), end = max(rpf$score)) # ?
one_iso_fractionLength <- fractionLength(one_iso_cds, sum(width(one_iso_exons)))
one_iso_disengagementScore <- disengagementScore(one_iso_cds, rpf, txdb)
one_iso_ribosomeReleaseScore <- ribosomeReleaseScore(one_iso_cds, rpf, threeUTRsByTranscript(txdb, use.names = TRUE), RNA = NULL)
one_iso_ribosomeStallingScore <- ribosomeStallingScore(one_iso_cds, rpf)
one_iso_orfScore <- orfScore(one_iso_cds, rpf)
#distToCds(one_iso, fiveUTRs, cds = NULL, extension = NULL)
one_iso_insideOutsideORF <- insideOutsideORF(one_iso_cds, rpf, txdb)

one_iso_cds_copy <- one_iso_cds
seqlevels(one_iso_cds_copy) <- gsub("^chr", "", seqlevels(one_iso_cds_copy))
ff <- FaFile("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
one_iso_kozakSequenceScore <- kozakSequenceScore(one_iso_cds_copy, ff, species = "human", include.N = FALSE)

one_iso_features <- data.table(tx = names(one_iso_cds),
                               fpkm = one_iso_fpkm,
                               entropy = one_iso_entropy,
                               floss = one_iso_floss,
                               fractionLength = one_iso_fractionLength,
                               disengagementScore = one_iso_disengagementScore,
                               ribosomeReleaseScore = one_iso_ribosomeReleaseScore,
                               ribosomeStallingScore = one_iso_ribosomeStallingScore,
                               orfScore = one_iso_orfScore$ORFScores,
                               insideOutsideORF = one_iso_insideOutsideORF,
                               kozakSequenceScore = one_iso_kozakSequenceScore)

one_iso_features$cdsLen <- sum(width(one_iso_cds))

### feature: % positions covered !!!!!!!!!!!!!!
oi_tiled <- tile1(one_iso_cds)

cov <- function(x, rpf) {
  co <- countOverlaps(x, rpf)
  cov <- sum(co > 0) / length(co)
  return(cov)
}

oi_coverage <- sapply(oi_tiled, function(x){cov(x, rpf)})

one_iso_features$coverage <- oi_coverage

# save(one_iso_features, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/isoform_deconvolution/one_iso_features.Rsave")





#####################################
### get off-frame ()

library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

seqlevels(genome)[which(seqlevels(genome) == "chrM")] <- "chrMT"

#gsub("^chrUn_", "", seqlevels(genome))
#gsub("v1$", ".1", seqlevels(genome))
#seqlevels(genome) <- gsub("^chrUn_", "", seqlevels(genome))
#seqlevels(genome) <- gsub("v1$", ".1", seqlevels(genome))

seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, one_iso_exons) # get whole tx

#minus <- seq[which(unlist(runValue(strand(one_iso_exons)) == "-"))]
#minus <- DNAStringSetList(sapply(minus, function(x){rev(x)}))
#seq[which(unlist(runValue(strand(one_iso_exons)) == "-"))] <- minus
seq <- DNAStringSet(sapply(seq, function(x){paste0(x, sep="", collapse="")}))
#seq_translated <- translate(seq)

# "ATG|TTG|CTG|AAG|ACG|AGG|ATA|ATT|ATC"
ORFs <- findMapORFs(one_iso_exons, seq, startCodon = "ATG|TTG|CTG",
            stopCodon = stopDefinition(1), longestORF = FALSE, minimumLength = 0)

ORFs_tx <- findORFs(seq, startCodon = "ATG|TTG|CTG",
                    stopCodon = stopDefinition(1), longestORF = FALSE, minimumLength = 0)
names(ORFs_tx) <- names(one_iso_exons)

# pick off-frame ORF(s) of similar length (check + and -)
one_iso_txlen <- txLengths[txLengths$tx_name %in% names(one_iso_exons),]
one_iso_txlen <- one_iso_txlen[order(one_iso_txlen$tx_name),] # get canonical CDS tx coordinates
# o_orfs <- o_orfs[(start(o_orfs) - can_start) %% 3 != 0] # off-frame

alt_orfs <- GRanges()
for (i in 1:length(one_iso_cds)) {
  gr_cds <- one_iso_cds[[i]]
  gr_exon <- one_iso_exons[[i]]
  gr_orfs_tx <- ORFs_tx[[i]]
  gr_orfs <- ORFs[[i]]
  can_orf_start <- one_iso_txlen[i,]$utr5_len + 1
  can_orf_end <- one_iso_txlen[i,]$utr5_len + one_iso_txlen[i,]$cds_len
  ## pick alternative ORFs
  # in-frame (extension/truncation)
  in_frame_idx <- which((start(gr_orfs_tx) - can_orf_start) %% 3 == 0)
  in_frame_idx <- in_frame_idx[sort(width(gr_orfs_tx[in_frame_idx]), index.return=TRUE, decreasing=TRUE)$ix[1:4]] # in-frame, top 4 longest
  # off-frame
  off_frame_idx <- which((start(gr_orfs_tx) - can_orf_start) %% 3 != 0)
  off_frame_idx <- off_frame_idx[sort(width(gr_orfs_tx[off_frame_idx]), index.return=TRUE, decreasing=TRUE)$ix[1:3]] # off-frame, top 3 longest
  # canonical CDS
  can_orf_idx <- which(start(gr_orfs_tx) == can_orf_start & end(gr_orfs_tx) == can_orf_end)
  # alternative ORFs indices, without canonical orf
  alt_idx <- c(in_frame_idx, off_frame_idx)
  if (length(which(in_frame_idx == can_orf_idx)) > 0) {
    alt_idx <- alt_idx[-which(in_frame_idx == can_orf_idx)]
  }
  
  ### retrieve alt orfs from ORFs (gr_orfs)
  for (a in alt_idx) {
    alt_orfs <- c(alt_orfs, gr_orfs[which(grepl(paste0("_", as.character(a), "$"), gr_orfs$names))])
  }
}

# split alt_orfs per names
alt_orfs <- split(alt_orfs, alt_orfs$names)

############## calculate features for alternative ORFs
alt_orfs_fpkm <- fpkm(alt_orfs, rpf)
alt_orfs_entropy <- entropy(alt_orfs, rpf)
alt_orfs_floss <- floss(alt_orfs, rpf, alt_orfs, start = min(rpf$score), end = max(rpf$score)) # ?
alt_orfs_fractionLength <- fractionLength(alt_orfs, rep(sum(width(one_iso_exons)), table(substr(names(alt_orfs), 1, 15))))
alt_orfs_disengagementScore <- disengagementScore(alt_orfs, rpf, txdb)
alt_orfs_ribosomeReleaseScore <- ribosomeReleaseScore(alt_orfs, rpf, threeUTRsByTranscript(txdb, use.names = TRUE), RNA = NULL)
alt_orfs_ribosomeStallingScore <- ribosomeStallingScore(alt_orfs, rpf)
alt_orfs_orfScore <- orfScore(alt_orfs, rpf)
alt_orfs_insideOutsideORF <- insideOutsideORF(alt_orfs, rpf, txdb)

alt_orfs_copy <- alt_orfs
seqlevels(alt_orfs_copy) <- gsub("^chr", "", seqlevels(alt_orfs_copy))
ff <- FaFile("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
alt_orfs_kozakSequenceScore <- kozakSequenceScore(alt_orfs_copy, ff, species = "human", include.N = FALSE)

alt_orfs_features <- data.table(tx = names(alt_orfs),
                               fpkm = alt_orfs_fpkm,
                               entropy = alt_orfs_entropy,
                               floss = alt_orfs_floss,
                               fractionLength = alt_orfs_fractionLength,
                               disengagementScore = alt_orfs_disengagementScore,
                               ribosomeReleaseScore = alt_orfs_ribosomeReleaseScore,
                               ribosomeStallingScore = alt_orfs_ribosomeStallingScore,
                               orfScore = alt_orfs_orfScore$ORFScores,
                               insideOutsideORF = alt_orfs_insideOutsideORF,
                               kozakSequenceScore = alt_orfs_kozakSequenceScore,
                               cdsLen = sum(width(alt_orfs)))


### feature: % positions covered !!!!!!!!!!!!!!
alt_orfs_tiled <- tile1(alt_orfs)
alt_orfs_tiled <- unlist(alt_orfs_tiled)
co <- countOverlaps(alt_orfs_tiled, rpf)
co <- split(co, names(co))

alt_orfs_coverage <- sapply(co, function(x){sum(x > 0) / length(x)})
alt_orfs_features$coverage <- alt_orfs_coverage

##### merge one_iso_features and alt_orf_features
one_iso_features$label <- rep("pos", nrow(one_iso_features))
alt_orfs_features$label <- rep("neg", nrow(alt_orfs_features))

one_iso_alt_orfs_features <- rbind(one_iso_features, alt_orfs_features)

# save(one_iso_alt_orfs_features, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/isoform_deconvolution/one_iso_alt_orfs_features.Rsave")

######### PLOT #########
ggplot(one_iso_alt_orfs_features, aes(x = fpkm, fill=label)) + geom_density(alpha=.3) + scale_x_log10() #+ theme(axis.title.x=element_text(size=40,face="bold"))
ggplot(one_iso_alt_orfs_features, aes(x = entropy, fill=label)) + geom_density(alpha=.3) + scale_x_log10()
ggplot(one_iso_alt_orfs_features, aes(x = floss, fill=label)) + geom_density(alpha=.3)
ggplot(one_iso_alt_orfs_features, aes(x = fractionLength, fill=label)) + geom_density(alpha=.3)
ggplot(one_iso_alt_orfs_features, aes(x = disengagementScore, fill=label)) + geom_density(alpha=.3) + scale_x_log10()
ggplot(one_iso_alt_orfs_features, aes(x = ribosomeReleaseScore, fill=label)) + geom_density(alpha=.3) + scale_x_log10()
ggplot(one_iso_alt_orfs_features, aes(x = ribosomeStallingScore, fill=label)) + geom_density(alpha=.3) + scale_x_log10()
ggplot(one_iso_alt_orfs_features, aes(x = orfScore, fill=label)) + geom_density(alpha=.3)
ggplot(one_iso_alt_orfs_features, aes(x = insideOutsideORF, fill=label)) + geom_density(alpha=.3) + scale_x_log10()
ggplot(one_iso_alt_orfs_features, aes(x = kozakSequenceScore, fill=label)) + geom_density(alpha=.3)

ggplot(one_iso_alt_orfs_features, aes(x=cdsLen, fill=label)) + geom_density(alpha=.3) + scale_x_log10()
ggplot(one_iso_alt_orfs_features, aes(x=coverage, fill=label)) + geom_density(alpha=.3)


# could pick extensions and truncations the same way
# but this way emphasizes coverage inside, periodicity - doesn't solve multiple isoform/overlap problem
# if I picked genes with 2 isoforms, each having unique part (either both or one with coverage)
# could try to make a model for decomposition based on them - modelling on one-isoform genes?








############################ TWO ISOFORMS ######################################
two_iso_cluster_numbers <- names(table(EnsemblClust$NA.)[table(EnsemblClust$NA.) == 2])
ensembl_two_iso <- EnsemblClust[EnsemblClust$NA. %in% two_iso_cluster_numbers]
seqlevels(ensembl_two_iso)[which(seqlevels(ensembl_two_iso) == "X")] <- "chrX"
seqlevels(ensembl_two_iso)[which(seqlevels(ensembl_two_iso) == "Y")] <- "chrY"
seqlevels(ensembl_two_iso)[which(seqlevels(ensembl_two_iso) == "MT")] <- "chrMT"
names(ensembl_two_iso) <- paste0(ensembl_two_iso$NA., ".", ensembl_two_iso$name)
# get only chromosomes (not scaffolds)
ensembl_two_iso <- ensembl_two_iso[grepl("^chr", as.character(seqnames(ensembl_two_iso)))]
ensembl_two_iso <- ensembl_two_iso[order(names(ensembl_two_iso))]
n <- ensembl_two_iso$name
cln <- names(ensembl_two_iso)
ensembl_two_iso <- blocks(ensembl_two_iso)
names(ensembl_two_iso) <- cln


# stumpf2013_SRR970538
ensembl_fpkm <- fpkm(ensembl_two_iso, rpf)
ensembl_fpkm <- data.table(tx = names(ensembl_two_iso), fpkm = ensembl_fpkm)
ntx <- unlist(strsplit(ensembl_fpkm$tx, split = "[.]"))
ncl <- ntx[seq(1, length(ntx), 2)]
ntx <- ntx[seq(2, length(ntx), 2)]
ensembl_fpkm$ncl <- ncl
ensembl_fpkm$ntx <- ntx
ensembl_fpkm <- ensembl_fpkm[order(ensembl_fpkm$ntx)]

# make sure they have 5'UTRs
tx_len <- txLengths[txLengths$tx_name %in% ensembl_fpkm$ntx,]
tx_len <- tx_len[order(tx_len$tx_name),]
ensembl_fpkm$utr5_len <- tx_len$utr5_len

ensembl_fpkm <- ensembl_fpkm[ensembl_fpkm$utr5_len > 0,]
ensembl_fpkm <- ensembl_fpkm[ensembl_fpkm$ncl %in% names(table(ensembl_fpkm$ncl)[table(ensembl_fpkm$ncl) == 2])] # 4190

two_iso_tx <- ensembl_fpkm[ensembl_fpkm$ncl %in% ensembl_fpkm[ensembl_fpkm$fpkm > 1]$ncl]$ntx # 1468 iso (734 genes)
two_iso_exons <- exons[names(exons) %in% two_iso_tx]
two_iso_exons <- two_iso_exons[order(names(two_iso_exons))] # 1468
two_iso_cds <- cds[names(cds) %in% two_iso_tx]
two_iso_cds <- two_iso_cds[order(names(two_iso_cds))] # 1468

seqlevels(two_iso_cds)[which(seqlevels(two_iso_cds) == "X")] <- "chrX"
seqlevels(two_iso_cds)[which(seqlevels(two_iso_cds) == "Y")] <- "chrY"
seqlevels(two_iso_cds)[which(seqlevels(two_iso_cds) == "MT")] <- "chrMT"
seqlevels(two_iso_exons)[which(seqlevels(two_iso_exons) == "X")] <- "chrX"
seqlevels(two_iso_exons)[which(seqlevels(two_iso_exons) == "Y")] <- "chrY"
seqlevels(two_iso_exons)[which(seqlevels(two_iso_exons) == "MT")] <- "chrMT"

# cluster in pairs, check unique exon parts - coverage on both or just one
names(two_iso_cds) <- paste0(ensembl_fpkm[ensembl_fpkm$ntx %in% names(two_iso_cds)]$ncl, ".", names(two_iso_cds))
names(two_iso_exons) <- paste0(ensembl_fpkm[ensembl_fpkm$ntx %in% names(two_iso_exons)]$ncl, ".", names(two_iso_exons))

two_iso_cds <- two_iso_cds[order(names(two_iso_cds))]
two_iso_exons <- two_iso_exons[order(names(two_iso_exons))]

# get every 2 (1st & 2nd; 3rd & 4th etc.)
a <- two_iso_cds[[99]]
b <- two_iso_cds[[100]]

#setdiff(a,b)
#intersect(a,b)

a <- GRanges(seqnames=Rle("chr1", 2), IRanges(start=c(1,7),end=c(4,12)), strand="+")
b <- GRanges(seqnames=Rle("chr1", 2), IRanges(start=c(3,10),end=c(6,15)), strand="+")

a_unique <- setdiff(a,b)
b_unique <- setdiff(b,a)
ab_intersect <- intersect(a,b)

# if either a or b are empty - test if the coverage over ab_intersect is significantly higher
countOverlaps(a_unique, rpf)
countOverlaps(b_unique, rpf)
countOverlaps(ab_intersect, rpf)

u <- 0
for (i in seq(1, length(two_iso_cds), 2)) {
  a <- two_iso_cds[[i]]
  b <- two_iso_cds[[i+1]]
  a_unique <- setdiff(a,b)
  b_unique <- setdiff(b,a)
  if (length(a_unique) > 0 & length(b_unique) > 0) {
    a_cov <- countOverlaps(a_unique, rpf)
    b_cov <- countOverlaps(b_unique, rpf)
    if (xor(a_cov == 0, b_cov == 0)) {
      if (a_cov == 0) {
        a_width <- sum(width(a_unique))
        if (a_width > 10) {
          u <- u + 1
          print(c("here is a:", i))
        }
      } else if (b_cov == 0) {
        b_width <- sum(width(b_unique))
        if (b_width > 10) {
          u <- u + 1
          print(c("here is b:", i))
        }
      }
    }
  }
}

print(u)


### if any has no-cov


iso_type <- c()
for (i in seq(1, length(two_iso_cds), 2)) {
  print(i)
  a <- two_iso_cds[[i]]
  b <- two_iso_cds[[i+1]]
  a_unique <- setdiff(a,b)
  b_unique <- setdiff(b,a)
  if (length(a_unique) == 0 & length(b_unique) == 0) { # exact
    iso_type <- c(iso_type, "exact_a", "exact_b")
  } else if (length(a_unique) > 0 & length(b_unique) > 0) { # both unique
    a_cov <- countOverlaps(a_unique, rpf)
    b_cov <- countOverlaps(b_unique, rpf)
    if (xor(a_cov == 0, b_cov == 0)) { # both unique, no cov on one
      if (a_cov == 0) {
        a_width <- sum(width(a_unique))
        if (a_width > 10) {
          iso_type <- c(iso_type, "both_unique_nocov_a", "both_unique_nocov_a")
        } else {
          iso_type <- c(iso_type, "both_unique_nocov_a_tooshort_a", "both_unique_nocov_a_tooshort_a")
        }
      } else if (b_cov == 0) {
        b_width <- sum(width(b_unique))
        if (b_width > 10) {
          iso_type <- c(iso_type, "both_unique_nocov_b", "both_unique_nocov_b")
        } else {
          iso_type <- c(iso_type, "both_unique_nocov_b_tooshort_b", "both_unique_nocov_b_tooshort_b")
        }
      }
    } else if (a_cov == 0 & b_cov == 0) { # both unique, no cov on both
      iso_type <- c(iso_type, "both_unique_nocov_both", "both_unique_nocov_both")
    } else if (a_cov > 0 & b_cov > 0) { # both unique, cov on both
      iso_type <- c(iso_type, "both_unique_cov_both", "both_unique_cov_both")
    }
  } else if (xor(length(a_unique) > 0, length(b_unique) > 0)) { # extension (one unique)
    if (length(a_unique) > 0) {
      a_cov <- countOverlaps(a_unique, rpf)
      if (a_cov == 0) { # extension a, no cov on a
        a_width <- sum(width(a_unique))
        if (a_width > 10) {
          iso_type <- c(iso_type, "extension_a_nocov_a", "extension_a_nocov_a")
        } else {
          iso_type <- c(iso_type, "extension_a_nocov_a_tooshort_a", "extension_a_nocov_a_tooshort_a")
        }
      } else { # extension a, cov on a
        iso_type <- c(iso_type, "extension_a_cov_a", "extension_a_cov_a")
      }
    } else if (length(b_unique) > 0) {
      b_cov <- countOverlaps(b_unique, rpf)
      if (b_cov == 0) { # extension b, no cov on b
        b_width <- sum(width(b_unique))
        if (b_width > 10) {
          iso_type <- c(iso_type, "extension_b_nocov_b", "extension_b_nocov_b")
        } else {
          iso_type <- c(iso_type, "extension_b_nocov_b_tooshort_b", "extension_b_nocov_b_tooshort_b")
        }
      } else { # extension b, cov on b
        iso_type <- c(iso_type, "extension_b_cov_b", "extension_b_cov_b")
      }
    }
  }
}

one_translated <- c("both_unique_nocov_a", "both_unique_nocov_b",
                    "extension_a_nocov_a", "extension_b_nocov_b")
decomposition <- c("both_unique_cov_both", "extension_a_cov_a", "extension_b_cov_b")

# save(one_iso_cds, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/isoform_deconvolution/one_iso_cds.Rsave")
# save(one_iso_exons, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/isoform_deconvolution/one_iso_exons.Rsave")
# save(two_iso_cds, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/isoform_deconvolution/two_iso_cds.Rsave")
# save(two_iso_exons, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/isoform_deconvolution/two_iso_exons.Rsave")
# save(iso_type, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/isoform_deconvolution/two_iso_type.Rsave")

### iso_type %in% one_translated (214: 107 genes)
### label +/-, calculate features (are any discriminating between the 2?)
### check intersection coverage comparing to unique (for decomposition)
it <- iso_type[iso_type %in% one_translated]
tic <- two_iso_cds[iso_type %in% one_translated]

labels <- c()
for (i in seq(1, length(tic), 2)) {
  desc <- it[i]
  a <- tic[[i]]
  b <- tic[[i+1]]
  if (desc == "both_unique_nocov_a") {
    labels <- c(labels, "neg", "pos")
  } else if (desc == "both_unique_nocov_b") {
    labels <- c(labels, "pos", "neg")
  } else if (desc == "extension_a_nocov_a") {
    labels <- c(labels, "neg", "pos")
  } else if (desc == "extension_b_nocov_b") {
    labels <- c(labels, "pos", "neg")
  }
}

two_iso_pos <- tic[labels == "pos"]
two_iso_neg <- tic[labels == "neg"]

np <- unlist(strsplit(names(two_iso_pos), split = "[.]"))
np <- np[seq(2, length(np), 2)]
names(two_iso_pos) <- np
two_iso_pos <- two_iso_pos[order(names(two_iso_pos))]

nn <- unlist(strsplit(names(two_iso_neg), split = "[.]"))
nn <- nn[seq(2, length(nn), 2)]
names(two_iso_neg) <- nn
two_iso_neg <- two_iso_neg[order(names(two_iso_neg))]

####### pos features
two_iso_pos_exons <- exons[names(two_iso_pos) %in% names(exons)]
two_iso_pos_exons <- two_iso_pos_exons[order(names(two_iso_pos_exons))]
### feature calculation
tip_fpkm <- fpkm(two_iso_pos, rpf)
tip_entropy <- entropy(two_iso_pos, rpf)
tip_floss <- floss(two_iso_pos, rpf, two_iso_pos, start = min(rpf$score), end = max(rpf$score)) # ?
tip_fractionLength <- fractionLength(two_iso_pos, sum(width(two_iso_pos_exons)))
tip_disengagementScore <- disengagementScore(two_iso_pos, rpf, txdb)
tip_ribosomeReleaseScore <- ribosomeReleaseScore(two_iso_pos, rpf, threeUTRsByTranscript(txdb, use.names = TRUE), RNA = NULL)
tip_ribosomeStallingScore <- ribosomeStallingScore(two_iso_pos, rpf)
tip_orfScore <- orfScore(two_iso_pos, rpf)
tip_insideOutsideORF <- insideOutsideORF(two_iso_pos, rpf, txdb)

two_iso_pos_copy <- two_iso_pos
seqlevels(two_iso_pos_copy) <- gsub("^chr", "", seqlevels(two_iso_pos_copy))
ff <- FaFile("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
tip_kozakSequenceScore <- kozakSequenceScore(two_iso_pos_copy, ff, species = "human", include.N = FALSE)

two_iso_pos_features <- data.table(tx = names(two_iso_pos),
                               fpkm = tip_fpkm,
                               entropy = tip_entropy,
                               floss = tip_floss,
                               fractionLength = tip_fractionLength,
                               disengagementScore = tip_disengagementScore,
                               ribosomeReleaseScore = tip_ribosomeReleaseScore,
                               ribosomeStallingScore = tip_ribosomeStallingScore,
                               orfScore = tip_orfScore$ORFScores,
                               insideOutsideORF = tip_insideOutsideORF,
                               kozakSequenceScore = tip_kozakSequenceScore,
                               label = rep("pos", length(two_iso_pos)),
                               cdsLen = sum(width(two_iso_pos)))


####### neg features
two_iso_neg_exons <- exons[names(two_iso_neg) %in% names(exons)]
two_iso_neg_exons <- two_iso_neg_exons[order(names(two_iso_neg_exons))]
### feature calculation
tin_fpkm <- fpkm(two_iso_neg, rpf)
tin_entropy <- entropy(two_iso_neg, rpf)
tin_floss <- floss(two_iso_neg, rpf, two_iso_neg, start = min(rpf$score), end = max(rpf$score)) # ?
tin_fractionLength <- fractionLength(two_iso_neg, sum(width(two_iso_neg_exons)))
tin_disengagementScore <- disengagementScore(two_iso_neg, rpf, txdb)
tin_ribosomeReleaseScore <- ribosomeReleaseScore(two_iso_neg, rpf, threeUTRsByTranscript(txdb, use.names = TRUE), RNA = NULL)
tin_ribosomeStallingScore <- ribosomeStallingScore(two_iso_neg, rpf)
tin_orfScore <- orfScore(two_iso_neg, rpf)
tin_insideOutsideORF <- insideOutsideORF(two_iso_neg, rpf, txdb)

two_iso_neg_copy <- two_iso_neg
seqlevels(two_iso_neg_copy) <- gsub("^chr", "", seqlevels(two_iso_neg_copy))
ff <- FaFile("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
tin_kozakSequenceScore <- kozakSequenceScore(two_iso_neg_copy, ff, species = "human", include.N = FALSE)

two_iso_neg_features <- data.table(tx = names(two_iso_neg),
                                   fpkm = tin_fpkm,
                                   entropy = tin_entropy,
                                   floss = tin_floss,
                                   fractionLength = tin_fractionLength,
                                   disengagementScore = tin_disengagementScore,
                                   ribosomeReleaseScore = tin_ribosomeReleaseScore,
                                   ribosomeStallingScore = tin_ribosomeStallingScore,
                                   orfScore = tin_orfScore$ORFScores,
                                   insideOutsideORF = tin_insideOutsideORF,
                                   kozakSequenceScore = tin_kozakSequenceScore,
                                   label = rep("neg", length(two_iso_neg)),
                                   cdsLen = sum(width(two_iso_neg)))

# merge feature tables
two_iso_features <- rbind(two_iso_pos_features, two_iso_neg_features)

ggplot(two_iso_features, aes(x = fpkm, fill=label)) + geom_density(alpha=.3) + scale_x_log10() #+ theme(axis.title.x=element_text(size=40,face="bold"))
ggplot(two_iso_features, aes(x = entropy, fill=label)) + geom_density(alpha=.3)
ggplot(two_iso_features, aes(x = floss, fill=label)) + geom_density(alpha=.3)
ggplot(two_iso_features, aes(x = fractionLength, fill=label)) + geom_density(alpha=.3)
ggplot(two_iso_features, aes(x = disengagementScore, fill=label)) + geom_density(alpha=.3) + scale_x_log10()
ggplot(two_iso_features, aes(x = ribosomeReleaseScore, fill=label)) + geom_density(alpha=.3) + scale_x_log10()
ggplot(two_iso_features, aes(x = ribosomeStallingScore, fill=label)) + geom_density(alpha=.3) + scale_x_log10()
ggplot(two_iso_features, aes(x = orfScore, fill=label)) + geom_density(alpha=.3)
ggplot(two_iso_features, aes(x = insideOutsideORF, fill=label)) + geom_density(alpha=.3) + scale_x_log10()
ggplot(two_iso_features, aes(x = kozakSequenceScore, fill=label)) + geom_density(alpha=.3)

ggplot(two_iso_features, aes(x=cdsLen, fill=label)) + geom_density(alpha=.3) + scale_x_log10()

ggplot(two_iso_features, aes(x = cdsLen, y = fpkm, colour=label)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(two_iso_features, aes(x = cdsLen, y = entropy, colour=label)) + geom_point(size=0.1) + scale_x_log10()
ggplot(two_iso_features, aes(x = cdsLen, y = floss, colour=label)) + geom_point(size=0.1) + scale_x_log10()
ggplot(two_iso_features, aes(x = cdsLen, y = fractionLength, colour=label)) + geom_point(size=0.1) + scale_x_log10()
ggplot(two_iso_features, aes(x = cdsLen, y = disengagementScore, colour=label)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(two_iso_features, aes(x = cdsLen, y = ribosomeReleaseScore, colour=label)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(two_iso_features, aes(x = cdsLen, y = ribosomeStallingScore, colour=label)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(two_iso_features, aes(x = cdsLen, y = orfScore, colour=label)) + geom_point(size=0.1) + scale_x_log10()
ggplot(two_iso_features, aes(x = cdsLen, y = insideOutsideORF, colour=label)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(two_iso_features, aes(x = cdsLen, y = kozakSequenceScore, colour=label)) + geom_point(size=0.1) + scale_x_log10()

### feature: % positions covered !!!!!!!!!!!!!!
tip_tiled <- tile1(two_iso_pos)

cov <- function(x, rpf) {
  co <- countOverlaps(x, rpf)
  cov <- sum(co > 0) / length(co)
  return(cov)
}

tip_coverage <- sapply(tip_tiled, function(x){cov(x, rpf)})

tin_tiled <- tile1(two_iso_neg)
tin_coverage <- sapply(tin_tiled, function(x){cov(x, rpf)})

two_iso_features$coverage <- c(tip_coverage, tin_coverage)

ggplot(two_iso_features, aes(x = coverage, fill=label)) + geom_density(alpha=.3)

# save(two_iso_features, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/isoform_deconvolution/two_iso_features.Rsave")

### iso_type %in% decomposition (932: 466 genes)
### try to decompose based on one iso intersection (check with RNA-seq level and riboMAP)

