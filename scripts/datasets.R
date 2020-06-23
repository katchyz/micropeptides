### datasets comparison, overlap of coding and non-coding
library(rtracklayer)
library(GenomicRanges)

ensembl <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/Ensembl/Homo_sapiens.GRCh38.79.chr.gtf.gz", format = "gtf")
tRNA <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/gencode/gencode.v28.tRNAs.gtf.gz", format = "gtf")

# check overlap with coding gene annotations
protein_coding <- c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "nonsense_mediated_decay", "non_stop_decay",
                    "polymorphic_pseudogene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "protein_coding")
pseudogene <- c("IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "processed_pseudogene",
                "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene",
                "translated_processed_pseudogene", "TR_J_pseudogene", "unprocessed_pseudogene")
long_noncoding <- c("3prime_overlapping_ncrna", "antisense", "lincRNA", "processed_transcript", "sense_intronic",
                    "sense_overlapping")
short_noncoding <- c("Mt_rRNA", "Mt_tRNA", "miRNA", "misc_RNA", "rRNA", "ribozyme", "sRNA", "scaRNA", "snRNA", "snoRNA")


ensembl_coding <- ensembl[ensembl$transcript_biotype %in% protein_coding]
ensembl_short_noncoding <- ensembl[ensembl$transcript_biotype %in% short_noncoding]

# remove non-coding annotations that overlap with coding regions
ensembl_short_noncoding <- ensembl_short_noncoding[-unique(to(findOverlaps(ensembl_coding, ensembl_short_noncoding)))]
tRNA <- tRNA[-unique(to(findOverlaps(ensembl_coding, tRNA)))]


# get riboseq
# ? get sequence, find ORFs








