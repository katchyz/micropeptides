### isoform deconvolution
#library(bedr) # cluster
library(rtracklayer)


# prior to loading: sort and cluster by bedtools (strand-specific)
# sort -k1,1 -k2,2n test.bed > test.sorted.bed
# bedtools cluster -s -i test.sorted.bed > test.clustered.bed

OpenProtClust <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/OpenProt/short_openprot.clustered.bed")
OpenProtClust <- unique(split(OpenProtClust, OpenProtClust$NA.)) # 14510
length(unlist(OpenProtClust)) # 15487

# joint Ensembl and RefSeq
EnsemblRefSeqClust <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/tx_ensembl_refseq.clustered.bed")
EnsemblRefSeqClust <- unique(split(EnsemblRefSeqClust, EnsemblRefSeqClust$NA.)) # 40236
length(unlist(EnsemblRefSeqClust)) # 201048

# Ensembl only
#EnsemblClust <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/Ensembl/tx_ensembl.clustered.bed")
#EnsemblClust <- unique(split(EnsemblClust, EnsemblClust$NA.)) # 20485

# RefSeq only
#RefSeqClust <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/RefSeq/tx_refseq.clustered.bed")
#RefSeqClust <- unique(split(RefSeqClust, RefSeqClust$NA.)) # 38268

# Fantom6
#fantom <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/fantom/F6_CAT.transcript.gtf.gz", format = "gtf")
#fantom <- split(fantom, fantom$transcript_id)
#export.bed(fantom, "/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/fantom/fantom.bed")
FantomClust <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/fantom/fantom.clustered.bed")
FantomClust <- FantomClust[order(FantomClust$name)]

fantom_info <- read.table(file="/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/fantom/fantom_IDs_type.tsv", sep="\t", header = TRUE)

l <- as.vector(sapply(as.character(fantom_info$trnscptID), function(x){length(unlist(strsplit(x, split=";")))}))
n <- as.character(unlist(sapply(as.character(fantom_info$trnscptID), function(x){unlist(strsplit(x, split=";"))})))

fantom_type <- data.table(type = rep(fantom_info$CAT_geneClass, l), tx = n)
fantom_type <- fantom_type[order(fantom_type$tx)]

FantomClust$type <- fantom_type$type

FantomClust <- unique(split(FantomClust, FantomClust$NA.)) # 75613
length(unlist(FantomClust)) # 631567

# leave clusters with at least one of: "lncRNA_antisense", "lncRNA_divergent", "lncRNA_intergenic", "lncRNA_sense_intronic"
# remove clusters with types: "coding_mRNA", "pseudogene", "uncertain_coding"
lnc_types <- c("lncRNA_antisense", "lncRNA_divergent", "lncRNA_intergenic", "lncRNA_sense_intronic")
coding_types <- c("coding_mRNA", "pseudogene", "uncertain_coding")

fantom_lnc <- list()
f <- 1
for (i in 1:length(FantomClust)) {
  type <- unique(FantomClust[[i]]$type)
  if (any(type %in% lnc_types)) {
    if (!any(type %in% coding_types)) {
      fantom_lnc[[f]] <- FantomClust[[i]]
      f <- f + 1
    }
  }
}

coding_idx <- sapply(FantomClust, function(x){any(unique(x$type) %in% coding_types)})
# remove coding, check for lncRNA
FantomLnc <- FantomClust[!coding_idx] # 48847
lnc_idx <- sapply(FantomLnc, function(x){any(unique(x$type) %in% lnc_types)})
FantomLnc <- FantomLnc[lnc_idx] # 44590

export.bed(unlist(FantomLnc), "/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/fantom/lnc_only.bed")
save(FantomLnc, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/fantom/lnc_only_clusters.Rsave")
# check ribo coverage

fan <- unlist(FantomLnc)
fn <- fan$name
fan <- blocks(fan)
names(fan) <- fn

bed_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/bed_tracks_for_ORFik"
beds <- list.files(path = bed_path)
i <- 95
rpf <- fread.bed(file.path(bed_path, beds[i]))
fan_fpkm <- fpkm(fan, rpf)
fan_fpkm <- data.table(tx = names(fan), fpkm = fan_fpkm)

#load(file = "/Volumes/USELESS/DATA/MICROPEPTIDES/OpenProt_short_gene_iso.Rsave") # isoforms - list (genes) of lists (isoforms)
#load(file = "/Volumes/USELESS/DATA/MICROPEPTIDES/OpenProt_short_isoforms.Rsave") # grl - flat list of isoforms
# export.bed(grl, "/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/OpenProt/openprot_short_gene_iso.bed")

# gene29: 8 isoforms
# gene357

If you just want to get the number of reads mapping to that junction
(not the actual reads), you can try featureCounts function in Rsubread
package. It has a parameter called 'countSplitAlignmentsOnly' which
allows you to count exon-spanning reads only. It is extremely fast.

library(Rsubread)
split_reads <- featureCounts(files = c("/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/human/Stern-Ginossar_N_2012.Human.fibroblasts.RPF.GRCh38.SRR592954.bam"),
              annot.ext = "/Volumes/USELESS/DATA/genomes/GTF/Homo_sapiens.GRCh38.79.chr.gtf",
              isGTFAnnotationFile = T, GTF.attrType = "transcript_id", splitOnly = T)





############ sORFs.org
sorfs_stumpf <- read.csv(file = "/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/sorfs_org/stumpf2013.txt", sep = "\t")


