################################################### TE
ribo_TE <- ribo_FPKM
ribo_TE <- ribo_TE[, -grep("^fritsch",colnames(ribo_TE)), with=FALSE]
ribo_TE <- ribo_TE[, -grep("^ingolia",colnames(ribo_TE)), with=FALSE]
ribo_TE <- ribo_TE[, -grep("^lee",colnames(ribo_TE)), with=FALSE]
ribo_TE <- ribo_TE[, -grep("^liu",colnames(ribo_TE)), with=FALSE]
ribo_TE <- ribo_TE[, -grep("^yoon",colnames(ribo_TE)), with=FALSE]
rownames(ribo_TE) <- ribo_FPKM$tx
ribo_TE$tx <- NULL

rna_TE <- rna_FPKM[, grep("^andreev",colnames(rna_FPKM)), with=FALSE]
rna_TE <- cbind(rna_TE, rna_FPKM[, grep("^gonzalez",colnames(rna_FPKM)), with=FALSE])
rna_TE <- cbind(rna_TE, rna_FPKM[,c("guo2010_SRR057513", "guo2010_SRR057514", "guo2010_SRR057518", "guo2010_SRR057519",
                                    "guo2010_SRR057523", "guo2010_SRR057524", "guo2010_SRR057527", "guo2010_SRR057530",
                                    "guo2010_SRR057533", "guo2010_SRR065776", "guo2010_SRR065777", "guo2010_SRR065781",
                                    "guo2010_SRR065782")])
rna_TE <- cbind(rna_TE, rna_FPKM[, grep("^hsieh",colnames(rna_FPKM)), with=FALSE])
rna_TE <- cbind(rna_TE, rna_FPKM[,c("rutkowski2015_SRR1523653", "rutkowski2015_SRR1523653", "rutkowski2015_SRR1523653",
                                    "rutkowski2015_SRR1523667", "rutkowski2015_SRR1523667", "rutkowski2015_SRR1523667")])
rna_TE <- cbind(rna_TE, rna_FPKM[, grep("^sidrauski",colnames(rna_FPKM)), with=FALSE])
rna_TE <- cbind(rna_TE, rna_FPKM[,c("stern2012_SRR592967", "stern2012_SRR592968", "stern2012_SRR592968",
                                    "stern2012_SRR592966", "stern2012_SRR592968", "stern2012_SRR592966",
                                    "stern2012_SRR592968", "stern2012_SRR592966", "stern2012_SRR592968",
                                    "stern2012_SRR592966")])
rna_TE <- cbind(rna_TE, rna_FPKM[, grep("^stumpf",colnames(rna_FPKM)), with=FALSE])
rna_TE <- cbind(rna_TE, rna_FPKM[,c("subtelny2014_SRR1039860")])
rownames(rna_TE) <- rna_FPKM$tx


ribo <- unlist(strsplit(colnames(ribo_TE), split = "_"))
ribo[seq(1, length(ribo), 2)]
ribo[seq(2, length(ribo), 2)]
rna <- unlist(strsplit(colnames(rna_TE), split = "_"))
rna[seq(2, length(rna), 2)]

matching_rna_ribo <- data.table(study = ribo[seq(1, length(ribo), 2)],
                                ribo = ribo[seq(2, length(ribo), 2)],
                                rna = rna[seq(2, length(rna), 2)])

save(matching_rna_ribo, file = "/Volumes/USELESS/DATA/MICROPEPTIDES/matching_rna_ribo.Rsave")

