# shift footprints in R
library(GenomicAlignments)
library(XML)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)

parseCigar <- function(cigar, shift, is_plus_strand) {
  c_signs <- unlist(explodeCigarOps(cigar))
  c_counts <- unlist(explodeCigarOpLengths(cigar))
  
  i = ifelse(is_plus_strand, 0, length(c_signs) + 1)
  increment = ifelse(is_plus_strand, 1, -1)
  limit = 0
  plusShift = 0
  while (shift >= limit) {
    i = i + increment
    if (c_signs[i] == "M") {
      limit = limit + c_counts[i]
    } else if (c_signs[i] == "N" || c_signs[i] == "D") {
      plusShift = plusShift + c_counts[i]
    } else if (c_signs[i] == "I") {
      plusShift = plusShift - c_counts[i]
    } else {
      warning(paste0("Not supported sign:", c_signs[i]))
    }
  }
  shift = shift + plusShift
  return(shift)
}

shiftFootprints <- function(footprints, selected_lengths, selected_shifts) {
  
  selected_shifts <- -1 * selected_shifts
  allFootrpintsShifted <- GRanges()
  
  if (length(selected_lengths) != length(selected_shifts)) {
    stop("Incorrect input. Not equal number of elements in",
         " selected_lengths and selected_shifts!")
  }
  if (sum(selected_lengths > abs(selected_shifts)) !=
      length(selected_shifts)) {
    stop("Incorrect input. selected_shifts cant be bigger",
         " than selected_lengths!")
  }
  
  for (i in 1:length(selected_lengths)) {
    message("Shifting footprints of length ", selected_lengths[i])
    riboReadsW <- footprints[qwidth(footprints) == selected_lengths[i]]
    if (length(riboReadsW) == 0) {
      next
    }
    is_cigar <- width(riboReadsW) != qwidth(riboReadsW)
    cigar_strings <- cigar(riboReadsW[is_cigar])
    sizes <- qwidth(riboReadsW)
    
    riboReadsW <- granges(riboReadsW, use.mcols = TRUE)
    riboReadsW$size <- sizes  #move footprint length to size
    riboReadsW <- resize(riboReadsW, 1)  #resize to 5' only
    
    cigars <- riboReadsW[is_cigar]
    notCigars <- riboReadsW[!is_cigar]
    
    # Not Cigars - shift 5' ends, + shift right, - shift left
    if (length(notCigars) != 0) {
      is_plus <- as.vector(strand(notCigars) == "+")
      shiftedNotCigarsPlus <- shift(notCigars[is_plus],
                                    selected_shifts[i])
      shiftedNotCigarsMinus <- shift(notCigars[!is_plus],
                                     -1 * selected_shifts[i])
      allFootrpintsShifted <- c(allFootrpintsShifted,
                                shiftedNotCigarsPlus,
                                shiftedNotCigarsMinus)
    }
    # Cigars
    if (length(cigars) != 0) {
      is_plus <- as.vector(strand(cigars) == "+")
      shift_by <- rep(selected_shifts[i], length(cigars))
      shift_by <- mapply(parseCigar, cigar_strings, shift_by, is_plus)
      shift_by[!is_plus] <- -1 * shift_by[!is_plus]
      shifted_cigars <- shift(cigars, shift_by)
      allFootrpintsShifted <- c(allFootrpintsShifted, shifted_cigars)
    }
  }
  
  message("Sorting shifted footprints...")
  allFootrpintsShifted <- sortSeqlevels(allFootrpintsShifted)
  allFootrpintsShifted <- sort(allFootrpintsShifted)
  return(allFootrpintsShifted)
}

# ## prepare empty granges
# txdb <- makeTxDbFromGFF("/export/valenfs/projects/fantom6/fancode/data/F6_CAT/F6_CAT.transcript.gtf", format = "gtf")
# txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
# # strand
# strand <- transcripts(txdb)
# strand <- data.frame(strand = strand(strand), tx_name = strand$tx_name)
# txLengths <- merge(txLengths, strand, by="tx_name", all=TRUE)
# rownames(txLengths) <- txLengths$tx_name
# txLengths <- txLengths[order(rownames(txLengths)),]


patterns <- c("Andreev.+",
              "Fritsch.+",
              "Gonzalez.+",
              "Guo.+",
              "Hsieh.+",
              "Ingolia_NT_2012.+",
              "Ingolia_NT_2014.+",
              "Lee.+",
              "Liu.+",
              "Rutkowski.+",
              "Sidrauski.+",
              "Stern.+",
              "Stumpf.+",
              "Subtelny.+",
              "Yoon.+")

offset_files <- list(c(rep("andreev1.xml", 2), rep("andreev2.xml", 4)),
              rep("fritsch.xml", 6),
              c("gonzalez1.xml", rep("gonzalez2.xml", 2), rep("gonzalez3.xml", 2)),
              rep("guo.xml" ,13),
              rep("hsieh.xml", 6),
              rep("ingolia2.xml", 3),
              rep("ingolia.xml", 4),
              c("lee3.xml", "lee2.xml", rep("lee1.xml", 3)),
              rep("liu.xml", 14),
              rep("rutkowski.xml", 6),
              rep("sidrauski.xml", 16),
              rep("stern.xml", 10),
              rep("stumpf.xml", 6),
              c("subtelny.xml"),
              rep("yoon.xml", 3))


save_prefixes <- c("andreev2015",
                   "fritsch2012",
                   "gonzalez2014",
                   "guo2010",  
                   "hsieh2012",    
                   "ingolia2012",  
                   "ingolia2014",
                   "lee2012",
                   "liu2013",    
                   "rutkowski2015",
                   "sidrauski2015",
                   "stern2012",
                   "stumpf2013",  
                   "subtelny2014",
                   "yoon2014")


####################################
# pattern <- "^Andreev.+"
# offsets_file <- "andreev1.xml"
####################################

offsets_path <- "/export/valenfs/projects/fantom6/fancode/data/F6_human_tracks/offsets"
bam_path = "/export/valenfs/projects/fantom6/fancode/data/F6_CAT_Ribo-seq_STARTalign"
libs <- list.files(path = bam_path, pattern = "\\.bam$")
shifted_granges_path <- "/export/valenfs/projects/fantom6/fancode/data/F6_human_tracks/tracks/"


for (p in 1:length(patterns)) {
  lp <- libs[grepl(patterns[p], libs)]
  srr <- strsplit(lp, "[.]")
  srr <- sapply(srr, function(x){strsplit(x[length(x)-3], "_")})
  srr <- sapply(srr, function(x){x[1]})
  for (l in 1:length(lp)) {
    # parse offsets
    x <- xmlParse(file.path(offsets_path, offset_files[[p]][l]))
    x <- data.frame(t(data.frame(xmlToList(x))))
    sel_lengths <- as.numeric(as.character(x$length))
    sel_shifts <- as.numeric(as.character(x$value))
    # load bam
    riboseq <- readGAlignments(file.path(bam_path, lp[l]))
    shifted <- shiftFootprints(riboseq, sel_lengths, sel_shifts)
    # save granges
    rdata_file <- paste0(shifted_granges_path, save_prefixes[p], "_", srr[l], ".RData")
    save(shifted, file = rdata_file)
    rm(riboseq)
    rm(shifted)
  }
}


#######################
# # parse offsets
# x <- xmlParse(file.path(offsets_path, offsets_file))
# x <- t(data.frame(xmlToList(x)))
# 
# sel_lengths <- as.numeric(x$length)
# sel_shifts <- as.numeric(x$value)
# 
# 
# libs <- list.files(path = wig_path, pattern = pattern)
# 
# riboseq <- readGAlignments("/Users/kasia/Documents/shoelaces/Data/example.bam")
# shiftFootprints(riboseq, c(28,29,30), c(-11,-12,-13))


