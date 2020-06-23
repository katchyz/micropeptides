## check gencode vs ensembl vs fantom lncRNA
library(rtracklayer)

gencode_basic <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/gencode.v28.basic.annotation.gtf.gz", format = "gtf")
gencode <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/gencode.v28.annotation.gtf.gz", format = "gtf")
ensembl <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/Homo_sapiens.GRCh38.79.chr.gtf.gz", format = "gtf")
fantom <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/F6_CAT.transcript.gtf.gz", format = "gtf")

fantom_info <- read.table(file="/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/fantom_IDs_type.tsv", sep="\t", header = TRUE)

tRNA <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/gencode.v28.tRNAs.gtf.gz", format = "gtf")

ccds <- import("/Volumes/USELESS/DATA/MICROPEPTIDES/annotations/CCDS.bed")
names(ccds) <- ccds$name
# ccds_len <- sapply(ccds$blocks, function(x){sum(width(x))}) # sum(ccds_len < 300) = 858

# gene
gb_gene <- unique(gencode_basic$gene_id)
gb_gene <- substr(gb_gene,1,15) # 58381
g_gene <- unique(gencode$gene_id)
g_gene <- substr(g_gene,1,15) # 58381
e_gene <- unique(ensembl$gene_id) # 65217
f_gene <- unique(fantom$gene_id) # 124047

sum(gb_gene %in% e_gene) # 56913
sum(g_gene %in% e_gene) # 56913
sum(e_gene %in% f_gene) # 37145

f_ensg <- f_gene[grepl(paste0("^", "ENSG"), f_gene)] # 38351
f_catg <- f_gene[grepl(paste0("^", "CATG"), f_gene)] # 85696

# tx
gb_tx <- unique(gencode_basic$transcript_id)
gb_tx <- substr(gb_tx,1,15) # 100986
g_tx <- unique(gencode$transcript_id)
g_tx <- substr(g_tx,1,15) # 203836
e_tx <- unique(ensembl$transcript_id) # 213623
f_tx <- unique(fantom$transcript_id) # 709176

sum(gb_tx %in% e_tx) # 96623
sum(g_tx %in% e_tx) # 194208
sum(e_tx %in% f_tx) # 0

f_ensg <- f_gene[grepl(paste0("^", "ENSG"), f_gene)] # 38351
f_catg <- f_gene[grepl(paste0("^", "CATG"), f_gene)] # 85696


### lincRNA

e_lincRNA <- unique(ensembl[ensembl$transcript_biotype == "lincRNA" & !is.na(ensembl$transcript_biotype),]$gene_id) # 8043


# sort(unique(ensembl$gene_biotype))
# [1] "3prime_overlapping_ncrna"           "IG_C_gene"                          "IG_C_pseudogene"                   
# [4] "IG_D_gene"                          "IG_J_gene"                          "IG_J_pseudogene"                   
# [7] "IG_V_gene"                          "IG_V_pseudogene"                    "Mt_rRNA"                           
# [10] "Mt_tRNA"                            "TEC"                                "TR_C_gene"                         
# [13] "TR_D_gene"                          "TR_J_gene"                          "TR_J_pseudogene"                   
# [16] "TR_V_gene"                          "TR_V_pseudogene"                    "antisense"                         
# [19] "lincRNA"                            "macro_lncRNA"                       "miRNA"                             
# [22] "misc_RNA"                           "non_coding"                         "polymorphic_pseudogene"            
# [25] "processed_pseudogene"               "processed_transcript"               "protein_coding"                    
# [28] "pseudogene"                         "rRNA"                               "ribozyme"                          
# [31] "sRNA"                               "scaRNA"                             "sense_intronic"                    
# [34] "sense_overlapping"                  "snRNA"                              "snoRNA"                            
# [37] "transcribed_processed_pseudogene"   "transcribed_unitary_pseudogene"     "transcribed_unprocessed_pseudogene"
# [40] "translated_processed_pseudogene"    "translated_unprocessed_pseudogene"  "unitary_pseudogene"                
# [43] "unprocessed_pseudogene"             "vaultRNA"


