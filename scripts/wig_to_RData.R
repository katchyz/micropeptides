### WIGs to granges (single)
# read single WIG, fwd, rev
# save to .RData
library(rtracklayer)

### import (separated by length?)
# andreev2015     "Andreev.+"
# fritsch2012     "Fritsch.+"
# gonzalez2014    "Gonzalez.+"
# guo2010         "Guo.+"
# hsieh2012       "Hsieh.+"
# ingolia2012     "Ingolia_NT_2012.+"
# ingolia2014     "Ingolia_NT_2014.+"
# lee2012         "Lee.+"
# liu2013         "Liu.+"
# rutkowski2015   "Rutkowski.+"
# sidrauski2015   "Sidrauski.+"
# stern2012       "Stern.+"
# stumpf2013      "Stumpf.+"
# subtelny2014    "Subtelny.+"
# yoon2014        "Yoon.+"

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


wig_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/tracks"

for (p in 1:length(patterns)) {
  pattern <- patterns[p]
  save_prefix <- save_prefixes[p]
  
  ##
  
  libs <- list.files(path = wig_path, pattern = pattern)
  srr <- unique(sapply(libs, function(x){substr(x,(nchar(x)-21),(nchar(x)-12))}))
  
  for (i in 1:length(srr)) {
    run = srr[i]
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
    ## save with run name
    rdata_file <- paste0("/Volumes/USELESS/DATA/MICROPEPTIDES/riboseq/", save_prefix, run, ".RData")
    save(wig, file = rdata_file)
  }
}


