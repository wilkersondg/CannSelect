# CannSelect probe design

library(stringr)
library(ggplot2)
library(reshape)
library(seqinr)
library(scales)

ThinBaitsByOverlap <- function(baits, ints) {
  FirstLast <- function(blist) {
      if (nrow(blist) == 1) {
        overlap <- NA
      } else {
        overlap <- blist$BAITEND[1] - blist$BAITSTART[nrow(blist)]
      }
      return(overlap)
  }
  KeepOrThin <- function(blist) {
    if (blist$INTLEN[1] > 1000) {
      blist$KORT <- "KEEP"
    } else {
      if (nrow(blist) == 1) {
        blist$KORT <- "KEEP"
      } else {
        if (nrow(blist) == 2) {
          if (blist$INTLEN[1] > 160 & blist$FLOVRLP[1] >= 100) {
            blist$KORT <- c("KEEP", "THIN")
          } else {
            blist$KORT <- c("KEEP", "KEEP")
          }
        } else {
          if (nrow(blist) == 3) {
            if (blist$FLOVRLP[1] > -70) {
              blist$KORT <- c("KEEP", "THIN", "KEEP")
            } else {
              blist$KORT <- c("KEEP", "KEEP", "KEEP")
            }
          } else {
            if (nrow(blist) == 4) {
              if (blist$FLOVRLP[1] > -70) {
                blist$KORT <- c("KEEP", "THIN", "THIN", "KEEP")
              } else {
                blist$KORT <- c("KEEP", "KEEP", "KEEP", "KEEP")
              }
            } else {
              if (nrow(blist) == 5) {
                if (blist$FLOVRLP[1] > -70) {
                  blist$KORT <- c("KEEP", "THIN", "THIN", "THIN", "KEEP")
                } else {
                  blist$KORT <- c("KEEP", "KEEP", "KEEP", "KEEP", "KEEP")
                }
              } else {
                if (nrow(blist) == 6) {
                  if ((round_any(abs(blist$FLOVRLP[1]), 140, f = round) / 140) == 1) {
                    if (sum(blist$BAITSTART[2:5] >= blist$BAITEND[1] &
                            blist$BAITEND[2:5] <= blist$BAITSTART[6]) == 1) {
                      kort <- c("THIN", "THIN", "THIN", "THIN")
                      kort[which(blist$BAITSTART[2:5] >= blist$BAITEND[1] &
                                   blist$BAITEND[2:5] <= blist$BAITSTART[6])] <- "KEEP"
                      blist$KORT <- c("KEEP", kort, "KEEP")
                    } else {
                      diffs <- c()
                      for (prb in 2:5) {
                        diffs <- append(diffs, abs(blist$BAITEND[1] - blist$BAITSTART[prb]) +
                          abs(blist$BAITEND[prb] - blist$BAITSTART[6]))
                      }
                      kort <- c("THIN", "THIN", "THIN", "THIN")
                      kort[which.min(diffs)] <- "KEEP"
                      blist$KORT <- c("KEEP", kort, "KEEP")
                    }
                  } else {
                      diffs1st <- c()
                      diffsLst <- c()
                      for (prb in 2:5) {
                        diffs1st <- append(diffs1st, abs(blist$BAITEND[1] - blist$BAITSTART[prb]))
                        diffsLst <- append(diffsLst, abs(blist$BAITEND[prb] - blist$BAITSTART[6]))
                      }
                      kort <- c("THIN", "THIN", "THIN", "THIN")
                      kort[unique(c(which.min(diffs1st), which.min(diffsLst)))] <- "KEEP"
                      blist$KORT <- c("KEEP", kort, "KEEP")
                  }
                } else {
                  if (nrow(blist) == 7) {
                    if ((round_any(abs(blist$FLOVRLP[1]), 140, f = round) / 140) == 1) {
                      if (sum(blist$BAITSTART[2:6] >= blist$BAITEND[1] &
                              blist$BAITEND[2:6] <= blist$BAITSTART[7]) == 1) {
                        kort <- c("THIN", "THIN", "THIN", "THIN", "THIN")
                        kort[which(blist$BAITSTART[2:6] >= blist$BAITEND[1] &
                                     blist$BAITEND[2:6] <= blist$BAITSTART[7])] <- "KEEP"
                        blist$KORT <- c("KEEP", kort, "KEEP")
                      } else {
                        diffs <- c()
                        for (prb in 2:6) {
                          diffs <- append(diffs, abs(blist$BAITEND[1] - blist$BAITSTART[prb]) +
                                            abs(blist$BAITEND[prb] - blist$BAITSTART[7]))
                        }
                        kort <- c("THIN", "THIN", "THIN", "THIN", "THIN")
                        kort[which.min(diffs)] <- "KEEP"
                        blist$KORT <- c("KEEP", kort, "KEEP")
                      }
                    } else {
                      diffs1st <- c()
                      diffsLst <- c()
                      for (prb in 2:6) {
                        diffs1st <- append(diffs1st, abs(blist$BAITEND[1] - blist$BAITSTART[prb]))
                        diffsLst <- append(diffsLst, abs(blist$BAITEND[prb] - blist$BAITSTART[7]))
                      }
                      kort <- c("THIN", "THIN", "THIN", "THIN", "THIN")
                      kort[unique(c(which.min(diffs1st), which.min(diffsLst)))] <- "KEEP"
                      blist$KORT <- c("KEEP", kort, "KEEP")
                    }
                  } else {
                    if (nrow(blist) == 8) {
                      if ((round_any(abs(blist$FLOVRLP[1]), 140, f = round) / 140) == 1) {
                        if (sum(blist$BAITSTART[2:7] >= blist$BAITEND[1] &
                                blist$BAITEND[2:7] <= blist$BAITSTART[8]) == 1) {
                          kort <- c("THIN", "THIN", "THIN", "THIN", "THIN", "THIN")
                          kort[which(blist$BAITSTART[2:7] >= blist$BAITEND[1] &
                                       blist$BAITEND[2:7] <= blist$BAITSTART[8])] <- "KEEP"
                          blist$KORT <- c("KEEP", kort, "KEEP")
                        } else {
                          diffs <- c()
                          for (prb in 2:7) {
                            diffs <- append(diffs, abs(blist$BAITEND[1] - blist$BAITSTART[prb]) +
                                              abs(blist$BAITEND[prb] - blist$BAITSTART[8]))
                          }
                          kort <- c("THIN", "THIN", "THIN", "THIN", "THIN", "THIN")
                          kort[which.min(diffs)] <- "KEEP"
                          blist$KORT <- c("KEEP", kort, "KEEP")
                        }
                      } else {
                        diffs1st <- c()
                        diffsLst <- c()
                        for (prb in 2:7) {
                          diffs1st <- append(diffs1st, abs(blist$BAITEND[1] - blist$BAITSTART[prb]))
                          diffsLst <- append(diffsLst, abs(blist$BAITEND[prb] - blist$BAITSTART[8]))
                        }
                        kort <- c("THIN", "THIN", "THIN", "THIN", "THIN", "THIN")
                        kort[unique(c(which.min(diffs1st), which.min(diffsLst)))] <- "KEEP"
                        blist$KORT <- c("KEEP", kort, "KEEP")
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    return(blist)
  }

  out <- NULL
  for (trg in seq(1:nrow(ints))) {
    prbnames <- unlist(str_split(ints$name[trg], "[|]"))
    baitnames <- unique(unlist(lapply(prbnames, function(x) {
      baits$TargetID[which(baits$Coordinates == x)]
    })))
    chr <- ifelse(ints$chrom[trg] == "NC_044371.1", "CHR01",
                        ifelse(ints$chrom[trg] == "NC_044375.1", "CHR02",
                               ifelse(ints$chrom[trg] == "NC_044372.1", "CHR03",
                                      ifelse(ints$chrom[trg] == "NC_044373.1", "CHR04",
                                             ifelse(ints$chrom[trg] == "NC_044374.1", "CHR05",
           ifelse(ints$chrom[trg] == "NC_044377.1", "CHR06",
                  ifelse(ints$chrom[trg] == "NC_044378.1", "CHR07",
                         ifelse(ints$chrom[trg] == "NC_044379.1", "CHR08",
                                ifelse(ints$chrom[trg] == "NC_044376.1", "CHR09",
                                       ifelse(ints$chrom[trg] == "NC_044370.1", "CHRX", "PROB"))))))))))
    numbait <- length(baitnames)
    blist <- as.data.frame(list("INTID" = rep(paste(chr, ":", ints$start[trg], "-", ints$end[trg], sep = ""), numbait),
                                "CHR" = rep(chr, numbait),
                                "INTSTART" = rep(ints$start[trg], numbait),
                                "INTEND" = rep(ints$end[trg], numbait),
                                "INTLEN" = rep(ints$length[trg], numbait),
                                "BAITNINT" = numbait,
                                "BAITID" = baitnames,
                                "BAITSTART" = as.numeric(word(word(baitnames, 2, sep = "[:]"), 1, sep = "[-]")),
                                "BAITEND" = as.numeric(word(word(baitnames, 2, sep = "[:]"), 2, sep = "[-]"))))
    blist <- blist[order(blist$BAITSTART), ]
    blist$FLOVRLP <- FirstLast(blist)
    blist$PNOVEL <- round((blist$INTLEN)/(blist$BAITNINT * 142)*100, 1)

    blist <- KeepOrThin(blist)

  out <- rbind(out, blist)
  }
return(out)
}

GetProbesFromSeq <- function(fastseq, meta) {

  out <- NULL

  for (seq in fastseq) {

    SEQID <- attr(seq, "name")

    pullmeta <- droplevels(meta[which(meta$transcript_accession == SEQID), ])

    annot <- unlist(str_split(getAnnot(seq)[[1]], " "))
    LOC <- annot[2]
    CHR <- str_replace(toupper(pullmeta$location), " ", "")
    CHR <- ifelse(CHR == "CHRX", CHR, str_replace(CHR, "CHR", "CHR0"))

    seq <- unlist(seq)

    cds <- substr(seq,
                  as.numeric(word(word(pullmeta$transcript_cds_coords, 2, sep = "[:]"), 1, sep = "[-]")),
                  as.numeric(word(word(pullmeta$transcript_cds_coords, 2, sep = "[:]"), 2, sep = "[-]")))

    TARGLEN <- str_count(cds)

    fsttile <- floor(TARGLEN/120)

    fstbeg <- c(seq.int(1, (fsttile*120), by = 120), (TARGLEN - 119))
    fstend <- c(seq.int(120, (fsttile*120), by = 120), TARGLEN)
    sndbeg <- c(seq.int(61, ((fsttile - 1)*120), by = 120), (TARGLEN - 179))
    sndend <- c(seq.int(180, (((fsttile - 1)*120) + 60), by = 120), (TARGLEN - 60))

    prbidx <- as.data.frame(list("START" = c(fstbeg, sndbeg), "END" = c(fstend, sndend)))
    prbidx <- prbidx[order(prbidx$START), ]

    prbidx$SEQ <- unlist(lapply(1:nrow(prbidx), function(x) {
      substr(seq, prbidx$START[x], prbidx$END[x])
    }))

    prbidx$SEQLEN <- str_count(prbidx$SEQ)

    prbidx$SEQID <- SEQID
    prbidx$TARGLEN <- TARGLEN
    prbidx$LOCID <- LOC
    prbidx$CHR <- CHR
    prbidx$REFCHR <- word(pullmeta$genomic_range, 1, sep = "[:]")

    prbidx$PRBID <- unlist(lapply(1:nrow(prbidx), function(x) {
      paste(prbidx$SEQID[1], "_CDS:", prbidx$START[x], "-", prbidx$END[x], sep = "")
    }))

    prbidx <- prbidx[ , c(9, 8, 7, 5, 6, 10, 1:4)]

    prbidx$PGC <- unlist(lapply(seq(1:nrow(prbidx)), function(x) {
      round((str_count(prbidx$SEQ[x], "G") +
               str_count(prbidx$SEQ[x], "C")) /
              str_length(prbidx$SEQ[x]) * 100, 2)
    }))

    prbidx$PREDREP <- ifelse(prbidx$PGC >= 44.56 & prbidx$PGC <= 66.25, "1",
                             ifelse((prbidx$PGC >= 42.92 & prbidx$PGC < 44.56) |
                                      (prbidx$PGC > 66.25 & prbidx$PGC <= 68.75), "2",
                                    ifelse((prbidx$PGC >= 40.42 & prbidx$PGC < 42.92) |
                                             (prbidx$PGC > 68.75 & prbidx$PGC <= 69.59), "3",
                                           ifelse((prbidx$PGC >= 38.75 & prbidx$PGC < 40.42) |
                                                    (prbidx$PGC > 69.59 & prbidx$PGC <= 70), "4",
                                                  ifelse(prbidx$PGC > 70, "4+",
                                                         ifelse(prbidx$PGC >= 34.59 & prbidx$PGC < 38.75, "5",
                                                                ifelse(prbidx$PGC >= 27.5 & prbidx$PGC < 34.59, "6",
                                                                       ifelse(prbidx$PGC < 27.5, "6+", "PROB"))))))))

    out <- rbind(out, prbidx)

  }
  return(out)
}

GetProbesFromMarkers <- function(markers, refgen) {

  out <- NULL

  for (chr in unique(marks$CHR)) {

    leftmark <- rightmark <- droplevels(marks[which(marks$CHR == chr), ])
    seqchr <- refgen[leftmark$REFCHR[1]]

    leftmark$START <- (leftmark$POS - 70)
    leftmark$END <- (leftmark$POS + 49)
    leftmark$PRBID <- paste(leftmark$REFCHR, "_", leftmark$POS, ":",
                            leftmark$START, "-", leftmark$END, sep = "")

    rightmark$START <- (rightmark$POS - 49)
    rightmark$END <- (rightmark$POS + 70)
    rightmark$PRBID <- paste(rightmark$REFCHR, "_", rightmark$POS, ":",
                             rightmark$START, "-", rightmark$END, sep = "")

    seqchr <- unlist(seqchr)

    markchr <- rbind(leftmark, rightmark)
    markchr <- markchr[order(markchr$START), ]

    markchr$SEQ <- unlist(lapply(1:nrow(markchr), function(x) {
      substr(seqchr, markchr$START[x], markchr$END[x])
    }))

    markchr$SEQLEN <- str_count(markchr$SEQ)

    markchr$PGC <- unlist(lapply(seq(1:nrow(markchr)), function(x) {
      round((str_count(markchr$SEQ[x], "G") + str_count(markchr$SEQ[x], "g") +
               str_count(markchr$SEQ[x], "C") + str_count(markchr$SEQ[x], "c")) /
              str_length(markchr$SEQ[x]) * 100, 2)
    }))

    markchr$PREDREP <- ifelse(markchr$PGC >= 44.56 & markchr$PGC <= 66.25, 1,
                              ifelse((markchr$PGC >= 42.92 & markchr$PGC < 44.56) |
                                       (markchr$PGC > 66.25 & markchr$PGC <= 68.75), 2,
                                     ifelse((markchr$PGC >= 40.42 & markchr$PGC < 42.92) |
                                              (markchr$PGC > 68.75 & markchr$PGC <= 69.59), 3,
                                            ifelse((markchr$PGC >= 38.75 & markchr$PGC < 40.42) |
                                                     (markchr$PGC > 69.59 & markchr$PGC <= 70), 4,
                                                   ifelse(markchr$PGC > 70, 4,
                                                          ifelse(markchr$PGC >= 34.59 & markchr$PGC < 38.75, 5,
                                                                 ifelse(markchr$PGC >= 27.5 & markchr$PGC < 34.59, 6,
                                                                        ifelse(markchr$PGC < 27.5, 6, NA))))))))

    markchr$TARGID <- word(markchr$PRBID, 1, sep = "[:]")

    markchr <- markchr[ , c(13, 1:12)]

    out <- rbind(out, markchr)

  }
  return(out)
}

Fill250KBGaps <- function(probes, markers) {

  probes <- droplevels(probes[which(probes$CHR != "CHRY"), ])
  gaplist <- NULL
  for (chr in unique(probes$CHR)) {
    pullchr <- droplevels(probes[which(probes$CHR == chr) , ])
    pullchr$UPSTR <- ifelse(pullchr$STRAND == "-", pullchr$REND, pullchr$RSTART)
    pullchr$DWNSTR <- ifelse(pullchr$STRAND == "-", pullchr$RSTART, pullchr$REND)
    pullchr <- pullchr[order(pullchr$UPSTR), ]
    pullchr$EOPP <- unlist(lapply(1:nrow(pullchr), function(prb) {
      if (prb == 1) {
        0
      } else {
        pullchr$DWNSTR[(prb - 1)]
      }
    }))
    pullchr$D2PP <- pullchr$UPSTR - pullchr$EOPP
    pullchr <- droplevels(pullchr[which(pullchr$D2PP >= 250000), ])
    pullchr$GAPID <- paste(chr, 1:nrow(pullchr), sep = "_")
    pullchr <- pullchr[ , c(10, 15, 13, 11, 14)]
    colnames(pullchr) <- c("CHR", "GAPID", "GAPSTART", "GAPEND", "GAPDIST")
    gaplist <- rbind(gaplist, pullchr)
  }

  gaplist$MARKNGAP <- unlist(lapply(seq(1:nrow(gaplist)), function(gap){
    pullchr <- droplevels(markers[which(markers$CHR == gaplist$CHR[gap]), ])
    pullchr <- pullchr[which(duplicated(pullchr$POS)), ]
    sum(pullchr$POS > gaplist$GAPSTART[gap] & pullchr$POS < gaplist$GAPEND[gap])
  }))
  gaplist$ROUND <- 0

  keptmarkers <- NULL

  for (chr in unique(gaplist$CHR)) {

    pullgap <- droplevels(gaplist[which(gaplist$CHR == chr), ])
    pullchr <- droplevels(markers[which(markers$CHR == chr), ])
    pullchr <- pullchr[which(duplicated(pullchr$POS)), ]

    keeppos <- NULL
    chrgap <- NULL

    for (gap in seq(1:nrow(pullgap))) {
      if (pullgap$MARKNGAP[gap] != 0) { # markers wn gap
        pullmark <- droplevels(pullchr[which(pullchr$POS > pullgap$GAPSTART[gap] &
                                               pullchr$POS < pullgap$GAPEND[gap]), ])
        first <- pullmark$POS[which.min(pullmark$POS)]
        last <- pullmark$POS[which.max(pullmark$POS)]
        dist <- last - first
        markngap <- sum(pullmark$POS > (first + 50000) &
                          pullmark$POS < (last - 50000))
        keeppos <- append(keeppos, c(first, last))
        newgap <- c(chr, paste(pullgap$GAPID[gap], ".1", sep = ""),
                    first, last, dist, markngap, 1)
        chrgap <- rbind(chrgap, newgap)
      }
    }

    chrgap <- as.data.frame(chrgap)
    colnames(chrgap) <- colnames(pullgap)
    for (c in 3:7) { chrgap[ , c] <- as.numeric(chrgap[ , c]) }
    gaplist <- rbind(gaplist, chrgap)

    whict <- 2

    while (sum(chrgap$GAPDIST > 250000 & chrgap$MARKNGAP != 0) > 0) {

      pullgap <- droplevels(chrgap[which(chrgap$GAPDIST > 250000 & chrgap$MARKNGAP != 0), ])
      chrgap <- NULL

      for (gap in seq(1:nrow(pullgap))) {
        pullmark <- droplevels(pullchr[which(pullchr$POS > (pullgap$GAPSTART[gap] + 50000) &
                                               pullchr$POS < (pullgap$GAPEND[gap] - 50000)), ])
        first <- pullmark$POS[which.min(pullmark$POS)]
        last <- pullmark$POS[which.max(pullmark$POS)]
        dist <- last - first
        markngap <- sum(pullmark$POS > (first + 50000) &
                          pullmark$POS < (last - 50000))
        keeppos <- append(keeppos, c(first, last))
        gid <- word(pullgap$GAPID[gap], 1, sep = "[.]")
        newgap <- c(chr, paste(gid, ".", whict, sep = ""),
                    first, last, dist, markngap, whict)
        chrgap <- rbind(chrgap, newgap)
      }

      chrgap <- as.data.frame(chrgap)
      colnames(chrgap) <- colnames(pullgap)
      for (c in 3:7) { chrgap[ , c] <- as.numeric(chrgap[ , c]) }
      gaplist <- rbind(gaplist, chrgap)

      whict <- whict + 1

    }

    keeppos <- keeppos[!duplicated(keeppos)]
    keeppos <- droplevels(markers[which(markers$CHR == chr & markers$POS %in% keeppos), ])

    keptmarkers <- rbind(keptmarkers, keeppos)

  }

  return(list(gaplist, keptmarkers))

}

###### SEQUENCE PROBES: GENE CDS TO FASTA FOR BLAST ######

cdsseq <- read.fasta("ncbi_dataset/data/rna.fna", seqtype = "DNA",
                     as.string = T, seqonly = F, forceDNAtolower = F)

cdsmeta <- read.table("ncbi_dataset/data/data_table.tsv", sep = "\t", header = T)[ , -c(4:6, 10, 12, 15:18)]
cdsmeta$gene_id <- as.factor(cdsmeta$gene_id) 

cdsmeta <- droplevels(cdsmeta[!grepl("XR", cdsmeta$transcript_accession), ])
cdsmeta <- droplevels(cdsmeta[!grepl("NM", cdsmeta$transcript_accession), ])

keepers <- unlist(lapply(unique(as.character(cdsmeta$gene_id)), function(gid) {
  pullgene <- cdsmeta[which(cdsmeta$gene_id == gid), ]
  if (nrow(pullgene) == 1) {
    pullgene$transcript_accession
  } else {
    pullgene$CDSLEN <- as.numeric(word(word(pullgene$transcript_cds_coords, 2, sep = "[:]"), 2, sep = "[-]")) -
      as.numeric(word(word(pullgene$transcript_cds_coords, 2, sep = "[:]"), 1, sep = "[-]"))
    if (length(pullgene$transcript_accession[which.max(pullgene$CDSLEN)]) > 1) {
      pullgene$transcript_accession[which.max(pullgene$CDSLEN)][1]
    } else {
      pullgene$transcript_accession[which.max(pullgene$CDSLEN)]
    }
  }
}))

cdsmeta <- droplevels(cdsmeta[which(cdsmeta$transcript_accession %in% keepers), ])
rm(keepers)

cdsmeta$CDSLEN <- as.numeric(word(word(cdsmeta$transcript_cds_coords, 2, sep = "[:]"), 2, sep = "[-]")) -
  as.numeric(word(word(cdsmeta$transcript_cds_coords, 2, sep = "[:]"), 1, sep = "[-]"))

cdsmeta$TWOXTILE <- (floor(cdsmeta$CDSLEN/120) + (floor(cdsmeta$CDSLEN/120) - 1) + 2)

cdsseq <- cdsseq[which(names(cdsseq) %in% cdsmeta$transcript_accession)]

cdsseq <- GetProbesFromSeq(cdsseq, cdsmeta)

write.table(cdsseq, "SeqProbeMeta.txt", quote = F, row.names = F, sep = "\t")

probeseqs <- as.list(cdsseq$SEQ)
names(probeseqs) <- cdsseq$PRBID

write.fasta(probeseqs, file.out = "SEQ.fasta", as.string = T,
            names = names(probeseqs), open = "w")

###### SEQUENCE PROBES: BLAST RESULTS ######

seqmeta <- read.table("SeqProbeMeta.txt", header = T)

blast <- read.table("SEQ_CS10v2.txt", header = T)

length(unique(blast$QSEQID))

blast$TARGID <- word(blast$QSEQID, 1, 2, sep = "_")
blast$JSCORE <- blast$LENGTH - (blast$MISMATCH + blast$GAPOPEN)

blast <- droplevels(blast[which(blast$JSCORE >= 60), ])

MM <- as.data.frame(table(blast$TARGID, blast$QSEQID))
MM <- droplevels(MM[which(MM$Freq > 0), ])

MM <- merge(MM, seqmeta[ , c(4, 6)],
            by.x = c("Var1", "Var2"), by.y = c("SEQID", "PRBID"), all = T)
MM$Freq <- ifelse(is.na(MM$Freq), 0, MM$Freq)

table(MM$Freq)

MM$KEEP <- ifelse((MM$Freq == 0 | MM$Freq > 3), 0, 1)

pkept <- aggregate(KEEP ~ Var1, data = MM, sum)
nprbs <- as.data.frame(table(MM$Var1))

pkept <- merge(pkept, nprbs, by = "Var1", all = T)
pkept$PKEPT <- round((pkept$KEEP / pkept$Freq) * 100, 2)

sum(pkept$PKEPT >= 90) 
sum(pkept$PKEPT >= 80) 
sum(pkept$PKEPT >= 70) 
sum(pkept$PKEPT >= 60) 
sum(pkept$PKEPT >= 50) 
sum(pkept$PKEPT >= 40) 
sum(pkept$PKEPT >= 30) 
sum(pkept$PKEPT >= 25) 

pkept <- merge(pkept, seqmeta[ , c(4:5)],
               by.x = "Var1", by.y = "SEQID", all.x = T, all.y = F)
pkept <- pkept[!duplicated(pkept), ]

par(mar = c(4.5, 4.5, 3, 1))
plot(x = pkept$PKEPT, y = pkept$TARGLEN,
     col = alpha("darkorchid4", .7), pch = 18, cex = 2,
     xlim = c(0, 102), ylim = c(0, 5650),
     xlab = "Coverage (%)", ylab = "Length (bp)",
     main = "Gene Target Length vs Coverage",
     xaxs = "i", yaxs = "i")

cdsmeta <- read.table("ncbi_dataset/data/data_table.tsv", sep = "\t", header = T)[ , -c(4:6, 10, 12, 15:18)]

pkept <- droplevels(pkept[which(pkept$PKEPT >= 70), ])

MM <- merge(MM, cdsmeta[ , c(7, 3)], by.x = "Var1", by.y = "transcript_accession", all.x = T, all.y = F)

MM_dropped <- droplevels(MM[which(!as.character(MM$Var1) %in% as.character(pkept$Var1)), ])
MM <- droplevels(MM[which(as.character(MM$Var1) %in% as.character(pkept$Var1)), ]) 
MM <- droplevels(MM[which(MM$KEEP == 1), ])

MM <- merge(MM, seqmeta[, c(6, 12)], by.x = "Var2", by.y = "PRBID", all.x = T, all.y = F)

seqmeta <- droplevels(seqmeta[which(seqmeta$PRBID %in% MM$Var2), ]) 

kept <- droplevels(blast[which(blast$QSEQID %in% MM$Var2), ]) 

bestblast <- lapply(unique(kept$QSEQID), function(x) {
  pullprb <- droplevels(kept[which(kept$QSEQID == x), ])
  if (nrow(pullprb) == 1) {
    as.character(pullprb)
  } else {
    as.character(pullprb[which.max(pullprb$BITSCORE), ])
  }
})
bestblast <- as.data.frame(do.call("rbind", bestblast)) 
colnames(bestblast) <- colnames(blast)
for (col in 3:12) {
  bestblast[ , col] <- as.numeric(bestblast[ , col])
}

sum(bestblast$PIDENT == 100 & bestblast$LENGTH == 120) 
sum(bestblast$PIDENT != 100 & bestblast$LENGTH == 120) 
sum(bestblast$PIDENT == 100 & bestblast$LENGTH != 120) 
sum(bestblast$PIDENT != 100 & bestblast$LENGTH != 120) 

modSE <- lapply(unique(bestblast$QSEQID), function(x) {
  pullprb <- droplevels(bestblast[which(bestblast$QSEQID == x), ])
  if (pullprb$LENGTH == 120) {
    c(pullprb$QSEQID, pullprb$SSEQID,
      pullprb$SSTART, pullprb$SEND)
  } else { # length != 120
    if ((pullprb$SEND - pullprb$SSTART) > 0) { # forward strand
      if (pullprb$QSTART == 1 & pullprb$QEND != 120) {
        c(pullprb$QSEQID, pullprb$SSEQID,
          pullprb$SSTART, (pullprb$SEND + (120 - pullprb$QEND)))
      } else {
        if (pullprb$QSTART != 1 & pullprb$QEND == 120) {
          c(pullprb$QSEQID, pullprb$SSEQID,
            ((pullprb$SSTART + 1) - (pullprb$QSTART)), pullprb$SEND)
        } else {
          if (pullprb$QSTART != 1 & pullprb$QEND != 120) {
            c(pullprb$QSEQID, pullprb$SSEQID,
              ((pullprb$SSTART + 1) - (pullprb$QSTART)), (pullprb$SEND + (120 - pullprb$QEND)))
          }
        }
      }
    } else { # reverse strand
      if (pullprb$QSTART == 1 & pullprb$QEND != 120) {
        c(pullprb$QSEQID, pullprb$SSEQID,
          pullprb$SSTART, (pullprb$SEND - (120 - pullprb$QEND)))
      } else {
        if (pullprb$QSTART != 1 & pullprb$QEND == 120) {
          c(pullprb$QSEQID, pullprb$SSEQID,
            ((pullprb$SSTART - 1) + (pullprb$QSTART)), pullprb$SEND)
        } else {
          if (pullprb$QSTART != 1 & pullprb$QEND != 120) {
            c(pullprb$QSEQID,  pullprb$SSEQID,
              ((pullprb$SSTART - 1) + (pullprb$QSTART)), (pullprb$SEND - (120 - pullprb$QEND)))
          }
        }
      }
    }
  }
})

modSE <- as.data.frame(do.call("rbind", modSE)) 
colnames(modSE) <- c("PRBID", "REFCHECK", "REFSTART", "REFEND")
for (col in 3:4) {
  modSE[ , col] <- as.numeric(modSE[ , col])
}

seqmeta <- merge(seqmeta, modSE, by = "PRBID", all = T) 

seqmeta <- droplevels(seqmeta[which(seqmeta$REFCHR == seqmeta$REFCHECK), ]) 

seqmeta <- merge(seqmeta, cdsmeta[ , c(3, 5, 7)],
                 by.x = "SEQID", by.y = "transcript_accession",
                 all.x = T, all.y = F)

write.table(seqmeta, "SeqPrbMeta_GRP2.txt", quote = F, row.names = F, sep = "\t")

### making gene list sup table ###
genetable <- merge(pkept, cdsmeta, by.x = "Var1", by.y = "transcript_accession", all.x = T, all.y = F)

write.table(genetable, "CannSelect_GeneList.txt", quote = F, row.names = F, sep = "\t")

###### MARKER PROBES: ADJUST HAPMAPS FOR MARKER MERGER ######

mckay <- read.table("MCKAY_CR85_POLY_140BP.hmp.txt", header = T, comment.char = "@")
mckay$rs. <- ifelse(mckay$chrom == 10, mckay$rs., str_replace(mckay$rs., "S", "S0"))
mckay$chrom <- as.character(ifelse(mckay$chrom == 10, mckay$chrom, paste("0", mckay$chrom, sep = "")))
write.table(mckay, "MCKAY_CR85_POLY_140BP_ADJ.hmp.txt", quote = F, row.names = F, sep = "\t")

gbs <- read.table("HEMPGBS_CR85_POLY_140BP.hmp.txt", header = T, comment.char = "@")
colnames(gbs)[24:ncol(gbs)] <- word(colnames(gbs)[24:ncol(gbs)], 1, sep = "[.]")
gbs$chrom <- ifelse(gbs$chrom == "NC_044371.1", "01",
                    ifelse(gbs$chrom == "NC_044375.1", "02",
                           ifelse(gbs$chrom == "NC_044372.1", "03",
                                  ifelse(gbs$chrom == "NC_044373.1", "04",
                                         ifelse(gbs$chrom == "NC_044374.1", "05",
                                                ifelse(gbs$chrom == "NC_044377.1", "06",
                                                       ifelse(gbs$chrom == "NC_044378.1", "07",
                                                              ifelse(gbs$chrom == "NC_044379.1", "08",
                                                                     ifelse(gbs$chrom == "NC_044376.1", "09",
                                                                            ifelse(gbs$chrom == "NC_044370.1", "10", "PROB"))))))))))
gbs <- droplevels(gbs[which(gbs$chrom != "PROB"), ])
gbs$rs. <- paste("S", gbs$chrom,"_", word(gbs$rs., 3, sep = "_"), sep = "")
write.table(gbs, "HEMPGBS_CR85_POLY_140BP_ADJ.hmp.txt", quote = F, row.names = F, sep = "\t")

woods <- read.table("WOODS_CR95_POLY_400BP.hmp.txt", header = T, comment.char = "@")
woods <- droplevels(woods[which(woods$chrom != "NC_029855.1"), ])
woods$chrom <- ifelse(woods$chrom == "NC_044371.1", "01",
                      ifelse(woods$chrom == "NC_044375.1", "02",
                             ifelse(woods$chrom == "NC_044372.1", "03",
                                    ifelse(woods$chrom == "NC_044373.1", "04",
                                           ifelse(woods$chrom == "NC_044374.1", "05",
                                                  ifelse(woods$chrom == "NC_044377.1", "06",
                                                         ifelse(woods$chrom == "NC_044378.1", "07",
                                                                ifelse(woods$chrom == "NC_044379.1", "08",
                                                                       ifelse(woods$chrom == "NC_044376.1", "09",
                                                                              ifelse(woods$chrom == "NC_044370.1", "10", "PROB"))))))))))
woods$rs. <- paste("S", woods$chrom,"_", woods$pos, sep = "")
write.table(woods, "WOODS_CR95_POLY_400BP_ADJ.hmp.txt", quote = F, row.names = F, sep = "\t")

ren <- read.table("REN21_CR100_POLY_140BP.hmp.txt", header = T, comment.char = "@")
ren$chrom <- ifelse(ren$chrom == 1, "10",
                      ifelse(ren$chrom == 2, "01",
                             ifelse(ren$chrom == 3, "03",
                                    ifelse(ren$chrom == 4, "04",
                                           ifelse(ren$chrom == 5, "05",
                                                  ifelse(ren$chrom == 6, "02",
                                                         ifelse(ren$chrom == 7, "09",
                                                                ifelse(ren$chrom == 8, "06",
                                                                       ifelse(ren$chrom == 9, "07",
                                                                              ifelse(ren$chrom == 10, "08", "PROB"))))))))))
ren$rs. <- paste("S", ren$chrom,"_", ren$pos, sep = "")
write.table(ren, "REN21_CR100_POLY_140BP_ADJ.hmp.txt", quote = F, row.names = F, sep = "\t")

###### MARKER PROBES: ADDING SOURCE TO MASTER MARKER LIST ######

ren <- read.table("REN21_CR100_POLY_140BP_ADJ.hmp.txt", header = T, comment.char = "@")[ , c(1, 3:4)]
mckay <- read.table("MCKAY_CR85_POLY_140BP_ADJ.hmp.txt", header = T, comment.char = "@")[ , c(1, 3:4)]
gbs <- read.table("HEMPGBS_CR85_POLY_140BP_ADJ.hmp.txt", header = T, comment.char = "@")[ , c(1, 3:4)]
woods <- read.table("WOODS_CR95_POLY_400BP_ADJ.hmp.txt", header = T, comment.char = "@")[ , c(1, 3:4)]
marks <- read.table("REN_MCKAY_GBS_WOODS_140BP.hmp.txt", header = T, comment.char = "@")[ , c(1, 3:4)]

colnames(ren) <- colnames(mckay) <- colnames(gbs) <- colnames(woods) <- colnames(marks) <- c("RS", "CHR", "POS")

marks$SOURCE <- ifelse((marks$RS %in% mckay$RS &
                          marks$RS %in% ren$RS &
                          marks$RS %in% gbs$RS &
                          marks$RS %in% woods$RS), "ALL",
                       ifelse((marks$RS %in% mckay$RS &
                                 marks$RS %in% ren$RS &
                                 (!marks$RS %in% gbs$RS) &
                                 (!marks$RS %in% woods$RS)), "MCKAY_REN",
                              ifelse((marks$RS %in% mckay$RS &
                                        (!marks$RS %in% ren$RS) &
                                        (marks$RS %in% gbs$RS) &
                                        (!marks$RS %in% woods$RS)), "MCKAY_GBS",
                                     ifelse((marks$RS %in% mckay$RS &
                                               (!marks$RS %in% ren$RS) &
                                               (!marks$RS %in% gbs$RS) &
                                               (marks$RS %in% woods$RS)), "MCKAY_WOODS",
                                            ifelse((!marks$RS %in% mckay$RS &
                                                      (marks$RS %in% ren$RS) &
                                                      (marks$RS %in% gbs$RS) &
                                                      (!marks$RS %in% woods$RS)), "REN_GBS",
                                                   ifelse((!marks$RS %in% mckay$RS &
                                                             (marks$RS %in% ren$RS) &
                                                             (!marks$RS %in% gbs$RS) &
                                                             (marks$RS %in% woods$RS)), "REN_WOODS",
                                                          ifelse((!marks$RS %in% mckay$RS &
                                                                    (!marks$RS %in% ren$RS) &
                                                                    (marks$RS %in% gbs$RS) &
                                                                    (marks$RS %in% woods$RS)), "GBS_WOODS",
                 ifelse((marks$RS %in% mckay$RS &
                           (!marks$RS %in% ren$RS) &
                           (!marks$RS %in% gbs$RS) &
                           (!marks$RS %in% woods$RS)), "MCKAY",
                        ifelse((!marks$RS %in% mckay$RS &
                                  (marks$RS %in% ren$RS) &
                                  (!marks$RS %in% gbs$RS) &
                                  (!marks$RS %in% woods$RS)), "REN",
                               ifelse((!marks$RS %in% mckay$RS &
                                         (!marks$RS %in% ren$RS) &
                                         (marks$RS %in% gbs$RS) &
                                         (!marks$RS %in% woods$RS)), "GBS",
                                      ifelse((!marks$RS %in% mckay$RS &
                                                (!marks$RS %in% ren$RS) &
                                                (!marks$RS %in% gbs$RS) &
                                                (marks$RS %in% woods$RS)), "WOODS",
                                             ifelse((marks$RS %in% mckay$RS &
                                                       (marks$RS %in% ren$RS) &
                                                       (marks$RS %in% gbs$RS) &
                                                       (!marks$RS %in% woods$RS)), "MCKAY_REN_GBS",
                                                    ifelse((marks$RS %in% mckay$RS &
                                                              (!marks$RS %in% ren$RS) &
                                                              (marks$RS %in% gbs$RS) &
                                                              (marks$RS %in% woods$RS)), "MCKAY_GBS_WOODS",
                                                           ifelse((marks$RS %in% mckay$RS &
                                                                     (marks$RS %in% ren$RS) &
                                                                     (!marks$RS %in% gbs$RS) &
                                                                     (marks$RS %in% woods$RS)), "MCKAY_REN_WOODS",
                                                                  ifelse((!marks$RS %in% mckay$RS &
                                                                            (marks$RS %in% ren$RS) &
                                                                            (marks$RS %in% gbs$RS) &
                                                                            (marks$RS %in% woods$RS)), "REN_GBS_WOODS", "PROB")))))))))))))))

marks$CHR <- ifelse(marks$CHR == 10, "CHRX", paste("CHR0", marks$CHR, sep = ""))

rm(gbs, mckay, ren, woods)

write.table(marks, "MARKS_MERGED.txt", quote = F, row.names = F, sep = "\t")

###### MARKER PROBES: FILTERING BASED ON PROXIMITY ######

marks <- read.table("MARKS_MERGED.txt", header = T) 

prbs <- read.table("PRBS_GP1_GP2SEQ.txt", header = T)

prbs$CHR <- ifelse(prbs$REFCHR == "NC_044371.1", "CHR01",
                   ifelse(prbs$REFCHR == "NC_044375.1", "CHR02",
                         ifelse(prbs$REFCHR == "NC_044372.1", "CHR03",
                                ifelse(prbs$REFCHR == "NC_044373.1", "CHR04",
                                       ifelse(prbs$REFCHR == "NC_044374.1", "CHR05",
                                              ifelse(prbs$REFCHR == "NC_044377.1", "CHR06",
                                                     ifelse(prbs$REFCHR == "NC_044378.1", "CHR07",
                                                            ifelse(prbs$REFCHR == "NC_044379.1", "CHR08",
                                                                   ifelse(prbs$REFCHR == "NC_044376.1", "CHR09",
                                                                          ifelse(prbs$REFCHR == "NC_044370.1", "CHRX", "PROB"))))))))))

prbs$CHR <- ifelse(is.na(prbs$CHR), "CHRY", prbs$CHR)

table(prbs$CHR)
# CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09  CHRX  CHRY

marks$START <- marks$POS - 70
marks$END <- marks$POS + 69

table(marks$CHR)
# CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09  CHRX

table(marks$SOURCE, marks$CHR)
#                 CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09  CHRX # TOTAL
# ALL                13     7     5     8     7     6     3     6     4     5 # 64
# GBS              2265  1410  1171  1680  1016  1148   864  1640  1225  1968 # 14387 *
# GBS_WOODS         182   134   106   128    80    94    70   156   102   146 # 1198
# MCKAY            8740 12907 12107 12196 11105 11032 10666  9688  9065 10717 # 108223 *
# MCKAY_GBS          88    67    76    74    56    41    40    93    54    62 # 651
# MCKAY_GBS_WOODS    57    33    41    38    31    24    14    46    22    19 # 325
# MCKAY_REN         570   901   826   857   712   668   579   689   599   673 # 7074
# MCKAY_REN_GBS      10    15    11    21     8     5     5    11     9    11 # 106
# MCKAY_REN_WOODS   180   146   171   179   134   115   117   177   126   118 # 1463
# MCKAY_WOODS      1357  1324  1192  1350  1098  1098   938  1229  1042  1073 # 11701
# REN             22740 30125 27506 28799 23542 25894 20816 24226 22441 24387 # 250476 *
# REN_GBS           177   116    75   135    70    67    79   133   115    92 # 1059
# REN_GBS_WOODS      35    26    16    19    20    18    18    32    20    12 # 216
# REN_WOODS        1387  1615  1455  1515  1288  1213  1115  1303  1100  1271 # 13262
# WOODS           13971 20050 18517 19900 16234 17584 14827 14969 13973 19498 # 169523 *

marks$OVRLAP <- unlist(lapply(1:nrow(marks), function(mkr) {
  chrfor <- droplevels(prbs[which(prbs$CHR == marks$CHR[mkr] & prbs$STRAND == "+"), ])
  chrrev <- droplevels(prbs[which(prbs$CHR == marks$CHR[mkr] & prbs$STRAND == "-"), ])
  if (sum(sum((marks$START[mkr] <= chrfor$RSTART & marks$END[mkr] >= chrfor$RSTART) |
              (marks$START[mkr] <= chrfor$REND & marks$END[mkr] >= chrfor$REND)),
          sum((marks$START[mkr] <= chrrev$REND & marks$END[mkr] >= chrrev$REND) |
              (marks$START[mkr] <= chrrev$RSTART & marks$END[mkr] >= chrrev$RSTART))) > 0) {
    "OVERLAP"
  } else {
    "NOPE"
  }
}))

table(marks$OVRLAP)
#   NOPE OVERLAP

marks$WN100KB <- unlist(lapply(1:nrow(marks), function(mkr) {
  if (marks$OVRLAP[mkr] == "OVERLAP") {
    "OVERLAP"
  } else {
    pullchr <- droplevels(prbs[which(prbs$CHR == marks$CHR[mkr]), ])
    prox <- abs(((pullchr$REND + pullchr$RSTART)/2) - marks$POS[mkr])
    if (prox[which.min(prox)] < 100000) {
      "WN100KB"
    } else {
      if (prox[which.min(prox)] < 250000) {
        "WN250KB"
      } else {
        if (prox[which.min(prox)] >= 250000 & prox[which.min(prox)] <= 500000) {
          "BN250&500KB"
        } else {
          "GRTR500KB"
        }
      }
    }
  }
}))

table(marks$WN100KB)
# BN250&500KB   GRTR500KB     OVERLAP     WN100KB     WN250KB

table(marks$CHR, marks$WN100KB)
#       BN250&500KB GRTR500KB OVERLAP WN100KB WN250KB
# CHR01        4028      4178     890   35260    7416
# CHR02        7670     12700     617   34971   12918
# CHR03       10762      4533     550   30689   16741
# CHR04        7113      9766     692   35241   14087
# CHR05        7684     11023     420   25187   11087
# CHR06        8596      5130     489   29868   14924
# CHR07        6630      7384     385   24357   11395
# CHR08        4324      4831     731   35015    9497
# CHR09        6546      5416     539   26547   10849
# CHRX         7651     14903     546   26396   10556

marks <- droplevels(marks[which(marks$WN100KB %in% c("BN250&500KB", "GRTR500KB")), ]) 

table(marks$SOURCE)
#           ALL             GBS       GBS_WOODS           MCKAY       MCKAY_GBS MCKAY_GBS_WOODS       MCKAY_REN
# MCKAY_REN_GBS MCKAY_REN_WOODS     MCKAY_WOODS             REN         REN_GBS   REN_GBS_WOODS       REN_WOODS
#         WOODS

write.table(marks, "Marks4Probes_GTR250KB.txt", quote = F, row.names = F, sep = "\t")

###### MARKER PROBES: PULLING MARKER PROBE REFERENCE SEQUENCES FOR BLAST ######

marks = read.table("Marks4Probes_GTR250KB.txt", header = T)[ , -c(5:8)]

marks$REFCHR <- ifelse(marks$CHR == "CHR01", "NC_044371.1",
                       ifelse(marks$CHR == "CHR02", "NC_044375.1",
                              ifelse(marks$CHR == "CHR03", "NC_044372.1",
                                     ifelse(marks$CHR == "CHR04", "NC_044373.1",
                                            ifelse(marks$CHR == "CHR05", "NC_044374.1",
                                                   ifelse(marks$CHR == "CHR06", "NC_044377.1",
                                                          ifelse(marks$CHR == "CHR07", "NC_044378.1",
                                                                 ifelse(marks$CHR == "CHR08", "NC_044379.1",
                                                                        ifelse(marks$CHR == "CHR09", "NC_044376.1",
                                                                               ifelse(marks$CHR == "CHRX", "NC_044370.1", "PROB"))))))))))

cs10 <- read.fasta("GCF_900626175.2_cs10_genomic.fasta", seqtype = "DNA",
                   as.string = T, seqonly = F, forceDNAtolower = F)
cs10 <- cs10[1:10]

markprbs <- GetProbesFromMarkers(marks, cs10) 

markprbs <- markprbs[which(markprbs$PGC > 20 & markprbs$PGC < 60), ] 
dups <- markprbs$TARGID[duplicated(markprbs$TARGID)]
markprbs <- markprbs[which(markprbs$TARGID %in% dups), ] 

markprbs <- markprbs[which(!grepl("N", markprbs$SEQ)), ]
dups <- markprbs$TARGID[duplicated(markprbs$TARGID)]
markprbs <- markprbs[which(markprbs$TARGID %in% dups), ]

markprbs$SEQ <- toupper(markprbs$SEQ)

table(markprbs[which(!duplicated(markprbs$TARGID)), ]$SOURCE)
#           ALL             GBS       GBS_WOODS           MCKAY       MCKAY_GBS MCKAY_GBS_WOODS
#     MCKAY_REN   MCKAY_REN_GBS MCKAY_REN_WOODS     MCKAY_WOODS             REN         REN_GBS
# REN_GBS_WOODS       REN_WOODS           WOODS

table(markprbs[which(!duplicated(markprbs$TARGID)), ]$CHR)
# CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09  CHRX

table(markprbs$CHR)
# CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09  CHRX

write.table(markprbs, "MarkerProbeMeta.txt", quote = F, row.names = F, sep = "\t")

# spliting by chr to reduce file size of blast results #

ch1 <- droplevels(markprbs[which(markprbs$CHR == "CHR01"), ])
ch2 <- droplevels(markprbs[which(markprbs$CHR == "CHR02"), ])
ch3 <- droplevels(markprbs[which(markprbs$CHR == "CHR03"), ])
ch4 <- droplevels(markprbs[which(markprbs$CHR == "CHR04"), ])
ch5 <- droplevels(markprbs[which(markprbs$CHR == "CHR05"), ])
ch6 <- droplevels(markprbs[which(markprbs$CHR == "CHR06"), ])
ch7 <- droplevels(markprbs[which(markprbs$CHR == "CHR07"), ])
ch8 <- droplevels(markprbs[which(markprbs$CHR == "CHR08"), ])
ch9 <- droplevels(markprbs[which(markprbs$CHR == "CHR09"), ])
chx <- droplevels(markprbs[which(markprbs$CHR == "CHRX"), ])

ch1seqs <- as.list(ch1$SEQ)
names(ch1seqs) <- ch1$PRBID

ch2seqs <- as.list(ch2$SEQ)
names(ch2seqs) <- ch2$PRBID

ch3seqs <- as.list(ch3$SEQ)
names(ch3seqs) <- ch3$PRBID

ch4seqs <- as.list(ch4$SEQ)
names(ch4seqs) <- ch4$PRBID

ch5seqs <- as.list(ch5$SEQ)
names(ch5seqs) <- ch5$PRBID

ch6seqs <- as.list(ch6$SEQ)
names(ch6seqs) <- ch6$PRBID

ch7seqs <- as.list(ch7$SEQ)
names(ch7seqs) <- ch7$PRBID

ch8seqs <- as.list(ch8$SEQ)
names(ch8seqs) <- ch8$PRBID

ch9seqs <- as.list(ch9$SEQ)
names(ch9seqs) <- ch9$PRBID

chxseqs <- as.list(chx$SEQ)
names(chxseqs) <- chx$PRBID

# CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09  CHRX

write.fasta(ch1seqs[1:8030], file.out = "MARK_CHR1-1.fasta", as.string = T, names = names(ch1seqs)[1:8030], open = "w")
write.fasta(ch1seqs[8031:16062], file.out = "MARK_CHR1-2.fasta", as.string = T, names = names(ch1seqs)[8031:16062], open = "w")

write.fasta(ch2seqs[1:13296], file.out = "MARK_CHR2-1.fasta", as.string = T, names = names(ch2seqs)[1:13296], open = "w")
write.fasta(ch2seqs[13297:26592], file.out = "MARK_CHR2-2.fasta", as.string = T, names = names(ch2seqs)[13297:26592], open = "w")
write.fasta(ch2seqs[26593:39888], file.out = "MARK_CHR2-3.fasta", as.string = T, names = names(ch2seqs)[26593:39888], open = "w")

write.fasta(ch3seqs[1:14972], file.out = "MARK_CHR3-1.fasta", as.string = T, names = names(ch3seqs)[1:14972], open = "w")
write.fasta(ch3seqs[14973:29944], file.out = "MARK_CHR3-2.fasta", as.string = T, names = names(ch3seqs)[14973:29944], open = "w")

write.fasta(ch4seqs[1:11034], file.out = "MARK_CHR4-1.fasta", as.string = T, names = names(ch4seqs)[1:11034], open = "w")
write.fasta(ch4seqs[11035:22066], file.out = "MARK_CHR4-2.fasta", as.string = T, names = names(ch4seqs)[11035:22066], open = "w")
write.fasta(ch4seqs[22067:33100], file.out = "MARK_CHR4-3.fasta", as.string = T, names = names(ch4seqs)[22067:33100], open = "w")

write.fasta(ch5seqs[1:12204], file.out = "MARK_CHR5-1.fasta", as.string = T, names = names(ch5seqs)[1:12204], open = "w")
write.fasta(ch5seqs[12205:24408], file.out = "MARK_CHR5-2.fasta", as.string = T, names = names(ch5seqs)[12205:24408], open = "w")
write.fasta(ch5seqs[24409:36612], file.out = "MARK_CHR5-3.fasta", as.string = T, names = names(ch5seqs)[24409:36612], open = "w")

write.fasta(ch6seqs[1:13416], file.out = "MARK_CHR6-1.fasta", as.string = T, names = names(ch6seqs)[1:13416], open = "w")
write.fasta(ch6seqs[13417:26834], file.out = "MARK_CHR6-2.fasta", as.string = T, names = names(ch6seqs)[13417:26834], open = "w")

write.fasta(ch7seqs[1:13712], file.out = "MARK_CHR7-1.fasta", as.string = T, names = names(ch7seqs)[1:13712], open = "w")
write.fasta(ch7seqs[13713:27424], file.out = "MARK_CHR7-2.fasta", as.string = T, names = names(ch7seqs)[13713:27424], open = "w")

write.fasta(ch8seqs[1:8950], file.out = "MARK_CHR8-1.fasta", as.string = T, names = names(ch8seqs)[1:8950], open = "w")
write.fasta(ch8seqs[8951:17902], file.out = "MARK_CHR8-2.fasta", as.string = T, names = names(ch8seqs)[8951:17902], open = "w")

write.fasta(ch9seqs[1:11730], file.out = "MARK_CHR9-1.fasta", as.string = T, names = names(ch9seqs)[1:11730], open = "w")
write.fasta(ch9seqs[11731:23462], file.out = "MARK_CHR9-2.fasta", as.string = T, names = names(ch9seqs)[11731:23462], open = "w")

write.fasta(chxseqs[1:11028], file.out = "MARK_CHRX-1.fasta", as.string = T, names = names(chxseqs)[1:11028], open = "w")
write.fasta(chxseqs[11029:22058], file.out = "MARK_CHRX-2.fasta", as.string = T, names = names(chxseqs)[11029:22058], open = "w")
write.fasta(chxseqs[22059:33089], file.out = "MARK_CHRX-3.fasta", as.string = T, names = names(chxseqs)[22059:33089], open = "w")
write.fasta(chxseqs[33090:44116], file.out = "MARK_CHRX-4.fasta", as.string = T, names = names(chxseqs)[33090:44116], open = "w")

###### MARKER PROBES: BLAST RESULTS BY CHR ######

#MARKERS
# CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09  CHRX

#PROBES
# CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09  CHRX

### CHR01 ###

dash1 <- read.table("MARK_CHR1-1_CS10v2.txt") # 2085965
dash2 <- read.table("MARK_CHR1-2_CS10v2.txt") # 3574928

blast <- rbind(dash1, dash2) # 5660893
rm(dash1, dash2)

colnames(blast) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH",
                     "MISMATCH", "GAPOPEN", "QSTART", "QEND",
                     "SSTART", "SEND", "EVALUE", "BITSCORE")

length(unique(blast$QSEQID)) # 16032

blast$JSCORE <- blast$LENGTH - (blast$MISMATCH + blast$GAPOPEN)
sum(blast$JSCORE < 60) # 263930

blast <- droplevels(blast[which(blast$JSCORE >= 60), ]) # 5396963

MM <- as.data.frame(table(blast$QSEQID)) # 16032
max(MM$Freq) # 14386

MM <- droplevels(MM[which(MM$Freq <= 3), ]) # 6850

MM$TARGID <- as.factor(word(MM$Var1, 1, sep = "[:]"))
dups <- as.character(MM$TARGID[duplicated(as.character(MM$TARGID))])
MM <- droplevels(MM[which(MM$TARGID %in% dups), ]) # P = 6554; M = 3277

MM$Freq <- NULL
colnames(MM)[1] <- "PRBID"
MM$CHR <- "CHR01"
MM <- MM[ , c(2, 1, 3)]

KEEP <- MM # 6554

rm(blast, dups, MM)

### CHR02 ###

dash1 <- read.table("MARK_CHR2-1_CS10v2.txt")
dash2 <- read.table("MARK_CHR2-2_CS10v2.txt")
dash3 <- read.table("MARK_CHR2-3_CS10v2.txt")

blast <- rbind(dash1, dash2, dash3) # 19502990
rm(dash1, dash2, dash3)

colnames(blast) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH",
                     "MISMATCH", "GAPOPEN", "QSTART", "QEND",
                     "SSTART", "SEND", "EVALUE", "BITSCORE")

length(unique(blast$QSEQID)) # 39835

blast$JSCORE <- blast$LENGTH - (blast$MISMATCH + blast$GAPOPEN)
sum(blast$JSCORE < 60) # 835718

blast <- droplevels(blast[which(blast$JSCORE >= 60), ]) # 18667272

MM <- as.data.frame(table(blast$QSEQID)) # 39835
max(MM$Freq) # 16818

MM <- droplevels(MM[which(MM$Freq <= 3), ]) # 14851

MM$TARGID <- as.factor(word(MM$Var1, 1, sep = "[:]"))
dups <- as.character(MM$TARGID[duplicated(as.character(MM$TARGID))])
MM <- droplevels(MM[which(MM$TARGID %in% dups), ]) # P = 13886; M = 6943

MM$Freq <- NULL
colnames(MM)[1] <- "PRBID"
MM$CHR <- "CHR02"
MM <- MM[ , c(2, 1, 3)]

KEEP <- rbind(KEEP, MM) # 20440

rm(blast, dups, MM)

### CHR03 ###

dash1 <- read.table("MARK_CHR3-1_CS10v2.txt")
dash2 <- read.table("MARK_CHR3-2_CS10v2.txt")

blast <- rbind(dash1, dash2) # 13706583
rm(dash1, dash2)

colnames(blast) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH",
                     "MISMATCH", "GAPOPEN", "QSTART", "QEND",
                     "SSTART", "SEND", "EVALUE", "BITSCORE")

length(unique(blast$QSEQID)) # 29882

blast$JSCORE <- blast$LENGTH - (blast$MISMATCH + blast$GAPOPEN)
sum(blast$JSCORE < 60) # 554572

blast <- droplevels(blast[which(blast$JSCORE >= 60), ]) # 13152011

MM <- as.data.frame(table(blast$QSEQID)) # 29882
max(MM$Freq) # 16827

MM <- droplevels(MM[which(MM$Freq <= 3), ]) # 12581
MM$TARGID <- as.factor(word(MM$Var1, 1, sep = "[:]"))
dups <- as.character(MM$TARGID[duplicated(as.character(MM$TARGID))])
MM <- droplevels(MM[which(MM$TARGID %in% dups), ]) # P = 11940; M = 5970

MM$Freq <- NULL
colnames(MM)[1] <- "PRBID"
MM$CHR <- "CHR03"
MM <- MM[ , c(2, 1, 3)]

KEEP <- rbind(KEEP, MM) # 32380

rm(blast, dups, MM)

### CHR04 ###

dash1 <- read.table("MARK_CHR4-1_CS10v2.txt")
dash2 <- read.table("MARK_CHR4-2_CS10v2.txt")
dash3 <- read.table("MARK_CHR4-3_CS10v2.txt")

blast <- rbind(dash1, dash2, dash3) # 21102535
rm(dash1, dash2, dash3)

colnames(blast) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH",
                     "MISMATCH", "GAPOPEN", "QSTART", "QEND",
                     "SSTART", "SEND", "EVALUE", "BITSCORE")

length(unique(blast$QSEQID)) # 33045

blast$JSCORE <- blast$LENGTH - (blast$MISMATCH + blast$GAPOPEN)
sum(blast$JSCORE < 60) # 1024322

blast <- droplevels(blast[which(blast$JSCORE >= 60), ]) # 20078213

MM <- as.data.frame(table(blast$QSEQID)) # 33045
max(MM$Freq) # 18619

MM <- droplevels(MM[which(MM$Freq <= 3), ]) # 11561
MM$TARGID <- as.factor(word(MM$Var1, 1, sep = "[:]"))
dups <- as.character(MM$TARGID[duplicated(as.character(MM$TARGID))])
MM <- droplevels(MM[which(MM$TARGID %in% dups), ]) # P = 10886; M = 5443

MM$Freq <- NULL
colnames(MM)[1] <- "PRBID"
MM$CHR <- "CHR04"
MM <- MM[ , c(2, 1, 3)]

KEEP <- rbind(KEEP, MM) # 43266

rm(blast, dups, MM)

### CHR05 ###

dash1 <- read.table("MARK_CHR5-1_CS10v2.txt")
dash2 <- read.table("MARK_CHR5-2_CS10v2.txt")
dash3 <- read.table("MARK_CHR5-3_CS10v2.txt")

blast <- rbind(dash1, dash2, dash3) # 17606789
rm(dash1, dash2, dash3)

colnames(blast) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH",
                     "MISMATCH", "GAPOPEN", "QSTART", "QEND",
                     "SSTART", "SEND", "EVALUE", "BITSCORE")

length(unique(blast$QSEQID)) # 36533

blast$JSCORE <- blast$LENGTH - (blast$MISMATCH + blast$GAPOPEN)
sum(blast$JSCORE < 60) # 670773

blast <- droplevels(blast[which(blast$JSCORE >= 60), ]) # 16936016

MM <- as.data.frame(table(blast$QSEQID)) # 36533
max(MM$Freq) # 17480

MM <- droplevels(MM[which(MM$Freq <= 3), ]) # 13322
MM$TARGID <- as.factor(word(MM$Var1, 1, sep = "[:]"))
dups <- as.character(MM$TARGID[duplicated(as.character(MM$TARGID))])
MM <- droplevels(MM[which(MM$TARGID %in% dups), ]) # P = 12470; M = 6235

MM$Freq <- NULL
colnames(MM)[1] <- "PRBID"
MM$CHR <- "CHR05"
MM <- MM[ , c(2, 1, 3)]

KEEP <- rbind(KEEP, MM) # 55736

rm(blast, dups, MM)

### CHR06 ###

dash1 <- read.table("MARK_CHR6-1_CS10v2.txt")
dash2 <- read.table("MARK_CHR6-2_CS10v2.txt")

blast <- rbind(dash1, dash2) # 13204123
rm(dash1, dash2)

colnames(blast) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH",
                     "MISMATCH", "GAPOPEN", "QSTART", "QEND",
                     "SSTART", "SEND", "EVALUE", "BITSCORE")

length(unique(blast$QSEQID)) # 26783

blast$JSCORE <- blast$LENGTH - (blast$MISMATCH + blast$GAPOPEN)
sum(blast$JSCORE < 60) # 510267

blast <- droplevels(blast[which(blast$JSCORE >= 60), ]) # 12693856

MM <- as.data.frame(table(blast$QSEQID)) # 26783
max(MM$Freq) # 16826

MM <- droplevels(MM[which(MM$Freq <= 3), ]) # 10298
MM$TARGID <- as.factor(word(MM$Var1, 1, sep = "[:]"))
dups <- as.character(MM$TARGID[duplicated(as.character(MM$TARGID))])
MM <- droplevels(MM[which(MM$TARGID %in% dups), ]) # P = 9730; M = 4865

MM$Freq <- NULL
colnames(MM)[1] <- "PRBID"
MM$CHR <- "CHR06"
MM <- MM[ , c(2, 1, 3)]

KEEP <- rbind(KEEP, MM) # 65446

rm(blast, dups, MM)

### CHR07 ###

dash1 <- read.table("MARK_CHR7-1_CS10v2.txt")
dash2 <- read.table("MARK_CHR7-2_CS10v2.txt")

blast <- rbind(dash1, dash2) # 11977303
rm(dash1, dash2)

colnames(blast) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH",
                     "MISMATCH", "GAPOPEN", "QSTART", "QEND",
                     "SSTART", "SEND", "EVALUE", "BITSCORE")

length(unique(blast$QSEQID)) # 27373

blast$JSCORE <- blast$LENGTH - (blast$MISMATCH + blast$GAPOPEN)
sum(blast$JSCORE < 60) # 518674

blast <- droplevels(blast[which(blast$JSCORE >= 60), ]) # 11458629

MM <- as.data.frame(table(blast$QSEQID))
max(MM$Freq) # 18236

MM <- droplevels(MM[which(MM$Freq <= 3), ]) # 9499
MM$TARGID <- as.factor(word(MM$Var1, 1, sep = "[:]"))
dups <- as.character(MM$TARGID[duplicated(as.character(MM$TARGID))])
MM <- droplevels(MM[which(MM$TARGID %in% dups), ]) # P = 8948; M = 4474

MM$Freq <- NULL
colnames(MM)[1] <- "PRBID"
MM$CHR <- "CHR07"
MM <- MM[ , c(2, 1, 3)]

KEEP <- rbind(KEEP, MM) # 74414

rm(blast, dups, MM)

### CHR08 ###

dash1 <- read.table("MARK_CHR8-1_CS10v2.txt")
dash2 <- read.table("MARK_CHR8-2_CS10v2.txt")

blast <- rbind(dash1, dash2) # 8160900
rm(dash1, dash2)

colnames(blast) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH",
                     "MISMATCH", "GAPOPEN", "QSTART", "QEND",
                     "SSTART", "SEND", "EVALUE", "BITSCORE")

length(unique(blast$QSEQID)) # 17877

blast$JSCORE <- blast$LENGTH - (blast$MISMATCH + blast$GAPOPEN)
sum(blast$JSCORE < 60) # 350045

blast <- droplevels(blast[which(blast$JSCORE >= 60), ]) # 7810855

MM <- as.data.frame(table(blast$QSEQID))
max(MM$Freq) # 15197

MM <- droplevels(MM[which(MM$Freq <= 3), ]) # 6268
MM$TARGID <- as.factor(word(MM$Var1, 1, sep = "[:]"))
dups <- as.character(MM$TARGID[duplicated(as.character(MM$TARGID))])
MM <- droplevels(MM[which(MM$TARGID %in% dups), ]) # P = 5946; M = 2973

MM$Freq <- NULL
colnames(MM)[1] <- "PRBID"
MM$CHR <- "CHR08"
MM <- MM[ , c(2, 1, 3)]

KEEP <- rbind(KEEP, MM) # 80360

rm(blast, dups, MM)

### CHR09 ###

dash1 <- read.table("MARK_CHR9-1_CS10v2.txt")
dash2 <- read.table("MARK_CHR9-2_CS10v2.txt")

blast <- rbind(dash1, dash2) #9703115
rm(dash1, dash2)

colnames(blast) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH",
                     "MISMATCH", "GAPOPEN", "QSTART", "QEND",
                     "SSTART", "SEND", "EVALUE", "BITSCORE")

length(unique(blast$QSEQID)) # 23404

blast$JSCORE <- blast$LENGTH - (blast$MISMATCH + blast$GAPOPEN)
sum(blast$JSCORE < 60) # 440725

blast <- droplevels(blast[which(blast$JSCORE >= 60), ]) # 9262390

MM <- as.data.frame(table(blast$QSEQID))
max(MM$Freq) # 13948

MM <- droplevels(MM[which(MM$Freq <= 3), ]) # 9715
MM$TARGID <- as.factor(word(MM$Var1, 1, sep = "[:]"))
dups <- as.character(MM$TARGID[duplicated(as.character(MM$TARGID))])
MM <- droplevels(MM[which(MM$TARGID %in% dups), ]) # P = 9148; M = 4574

MM$Freq <- NULL
colnames(MM)[1] <- "PRBID"
MM$CHR <- "CHR09"
MM <- MM[ , c(2, 1, 3)]

KEEP <- rbind(KEEP, MM) # 89508

rm(blast, dups, MM)

### CHRX ###

dash1 <- read.table("MARK_CHRX-1_CS10v2.txt")
dash2 <- read.table("MARK_CHRX-2_CS10v2.txt")
dash3 <- read.table("MARK_CHRX-3_CS10v2.txt")
dash4 <- read.table("MARK_CHRX-4_CS10v2.txt")

blast <- rbind(dash1, dash2, dash3, dash4) # 19722145
rm(dash1, dash2, dash3, dash4)

colnames(blast) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH",
                     "MISMATCH", "GAPOPEN", "QSTART", "QEND",
                     "SSTART", "SEND", "EVALUE", "BITSCORE")

length(unique(blast$QSEQID)) # 44030

blast$JSCORE <- blast$LENGTH - (blast$MISMATCH + blast$GAPOPEN)
sum(blast$JSCORE < 60) # 671547

blast <- droplevels(blast[which(blast$JSCORE >= 60), ]) # 19050598

MM <- as.data.frame(table(blast$QSEQID))
max(MM$Freq) # 16717

MM <- droplevels(MM[which(MM$Freq <= 3), ]) # 15637
MM$TARGID <- as.factor(word(MM$Var1, 1, sep = "[:]"))
dups <- as.character(MM$TARGID[duplicated(as.character(MM$TARGID))])
MM <- droplevels(MM[which(MM$TARGID %in% dups), ]) # P = 14806; M = 7403

MM$Freq <- NULL
colnames(MM)[1] <- "PRBID"
MM$CHR <- "CHRX"
MM <- MM[ , c(2, 1, 3)]

KEEP <- rbind(KEEP, MM) # 104314

rm(blast, dups, MM)

table(KEEP$CHR)
# CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09  CHRX

### SUBSET MARK META FOR GAPFILL STEP ###

markmeta <- read.table("MarkerProbeMeta.txt", header = T) # 295344

keepmeta <- droplevels(markmeta[which(markmeta$PRBID %in% KEEP$PRBID), ]) # 104314

write.table(keepmeta, "MkrPrbMeta_BLST.txt", quote = F, row.names = F, sep = "\t")

table(keepmeta$CHR)
# CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09  CHRX

table(keepmeta$SOURCE)
table(keepmeta$SOURCE, keepmeta$CHR)
#                 CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09 CHRX TOTAL
# ALL                 0     2     0     0     0     0     0     0     0    0     2
# GBS               284   290   230   292   184   148   150   156   174 1052  2960
# GBS_WOODS          18    32    34    16    18    24    20    30    16   58   266
# MCKAY            1250  2876  2502  2302  2658  1748  2092  1192  1946 2258 20824
# MCKAY_GBS           2     0     4     2     0     0     2     0     0    0    10
# MCKAY_GBS_WOODS     2     2     4     0     4     0     0     0     2    2    16
# MCKAY_REN          66   188   200   140   210   116   102    72   150  130  1374
# MCKAY_REN_GBS       0     0     2     2     0     0     0     0     0    0     4
# MCKAY_REN_WOODS    20    32    42    36    34    20    30    20    24   24   282
# MCKAY_WOODS       164   352   290   266   262   186   178   124   204  150  2176
# REN              2838  6262  5478  4872  5424  4816  3928  2642  4302 6478 47040
# REN_GBS             8    24    16    36    12     8    16    14    26   44   204
# REN_GBS_WOODS       0     6     6     2     0     4     6     2     4    4    34
# REN_WOODS         164   382   336   260   332   194   222   168   240  294  2592
# WOODS            1738  3438  2796  2660  3332  2466  2202  1526  2060 4312 26530

###### MARKER PROBES: SELECT MARKERS TO FILL IN GAPS ######

markmeta <- read.table("MkrPrbMeta_BLST.txt", header = T) # 104314

prbs <- read.table("PRBS_GP1_GP2SEQ.txt", header = T)

prbs$CHR <- ifelse(prbs$REFCHR == "NC_044371.1", "CHR01",
                   ifelse(prbs$REFCHR == "NC_044375.1", "CHR02",
                          ifelse(prbs$REFCHR == "NC_044372.1", "CHR03",
                                 ifelse(prbs$REFCHR == "NC_044373.1", "CHR04",
                                        ifelse(prbs$REFCHR == "NC_044374.1", "CHR05",
                                               ifelse(prbs$REFCHR == "NC_044377.1", "CHR06",
                                                      ifelse(prbs$REFCHR == "NC_044378.1", "CHR07",
                                                             ifelse(prbs$REFCHR == "NC_044379.1", "CHR08",
                                                                    ifelse(prbs$REFCHR == "NC_044376.1", "CHR09",
                                                                           ifelse(prbs$REFCHR == "NC_044370.1", "CHRX", "PROB"))))))))))

prbs$CHR <- ifelse(is.na(prbs$CHR), "CHRY", prbs$CHR)

table(prbs$CHR)
# CHR01 CHR02 CHR03 CHR04 CHR05 CHR06 CHR07 CHR08 CHR09  CHRX  CHRY

catch <- Fill250KBGaps(prbs, markmeta)

gapmeta <- catch[[1]]
grp2mkr <- catch[[2]]

gapmeta$GAP <- word(gapmeta$GAPID, 1, sep = "[.]")

ggplot(gapmeta, aes(x = as.factor(ROUND), y = GAPDIST)) +
  geom_boxplot()

lastround <- unlist(lapply(unique(gapmeta$GAP), function(gap) {
  pullgap <- droplevels(gapmeta[which(gapmeta$GAP == gap), ])
  pullgap$GAPID[which.max(pullgap$ROUND)]
  }))

gapmeta$FINAL <- ifelse(gapmeta$GAPID %in% lastround, "FINAL", "NOPE")

gapmeta$GAP <- as.factor(gapmeta$GAP)

write.table(gapmeta, "GP2MKR_GAPMETA.txt", quote = F, row.names = F, sep = "\t")
write.table(grp2mkr, "GP2MKR_PROBES.txt", quote = F, row.names = F, sep = "\t")

grp2mkr$GRP <- "GRP2SNP"
grp2mkr$STRAND <- "+"

grp2mkr <- grp2mkr[ , c(1, 9, 10, 13, 14, 6, 7, 8, 15, 3)]
colnames(grp2mkr) <- c("TARGID", "PRBID", "SEQ", "REP",
                       "GRP", "REFCHR", "RSTART", "REND",
                       "STRAND", "CHR")

prbs <- rbind(prbs, grp2mkr) # 18237

prbs$GRP <- ifelse(prbs$GRP %in% c("CHRY", "SYNSEQ"), "GRP1SEQ", prbs$GRP)

table(prbs$GRP)
# GRP1SEQ GRP1SNP GRP2SEQ GRP2SNP

sum(prbs$REP) # 49877

prbs <- prbs[order(prbs$CHR, prbs$RSTART), ]

write.table(prbs, "ALLPROBES.txt", quote = F, row.names = F, sep = "\t")

###### AGILENT BLAT SPECIFICITY CHECK ######

probes <- read.table("ALLPROBES.txt", header = T) # 18237
probes$PRBID <- str_replace(probes$PRBID, "Probe_", "")
probes$PRBID <- paste("PID", probes$PRBID, sep = "_")

hiblat <- read.table("BLAThitsOver5_CornellHempV3.txt", header = T) # 213
hiblat <- droplevels(hiblat[grepl("PID_NC", hiblat$ID), ]) # 150
hiblat$SOURCE <- unlist(lapply(hiblat$ID, function(x) {
  probes$GRP[which(probes$PRBID == x)]
  }))

grp2 <- probes[which(probes$GRP == "GRP2SNP"), ] # 5468
nogrp2 <- probes[which(probes$GRP != "GRP2SNP"), ] # 12769

grp2 <- droplevels(grp2[which(!(grp2$PRBID %in% hiblat$ID)), ]) # 5318
dups <- as.character(grp2$TARGID[duplicated(as.character(grp2$TARGID))])
grp2 <- droplevels(grp2[which(grp2$TARGID %in% dups), ]) # 5275

probes <- rbind(nogrp2, grp2) # 18043
probes <- probes[order(probes$CHR, probes$REFCHR), ]

write.table(probes, "ALLPROBES.txt", quote = F, row.names = F, sep = "\t")

# Final four-column SureDesign export
###### 4-COLUMN FORMAT FOR SUREDESIGN ######

probes <- read.table("ALLPROBES.txt", header = T) # 18043

probes <- probes[ , 1:4]
probes$REP <- 1

probes$PRBID <- str_replace(probes$PRBID, "Probe_", "")
probes$PRBID <- str_replace(probes$PRBID, "ProbeID_", "")
#probes$PRBID <- paste("PID", probes$PRBID, sep = "_")

probes <- droplevels(probes[!grepl("N", probes$SEQ), ]) # 18040

colnames(probes) <- c("TargetID", "ProbeID", "Sequence", "Replication")

write.table(probes, "CornellHempRedesign-4Column.txt", quote = F, row.names = F, sep = "\t")
# Final probe interval list
###### CREATING INTERVAL LIST FOR HSMETRICS ######

probes <- read.table("ALLPROBES.txt", header = T) # 18043
final <- read.table("CornellHempDesign1v3_Probes.txt", header = T, sep = "\t") # 18040

length(final$ProbeID[which(!(final$ProbeID %in% probes$PRBID))]) # 105
# These are the probes designed around the CBD Synthases; list makes up the discrepancy below

probes <- droplevels(probes[which(probes$PRBID %in% final$ProbeID), ]) # 17935

probes <- droplevels(probes[which(!is.na(probes$REFCHR)), ]) # 17935

probes$STRAND <- ifelse(probes$REND < probes$RSTART, "-", "+")
sum(probes$STRAND == "-") # 1213

# format is weird, coordinates have to be in ascending order
# start must be 'lower' than end; regardless of strand
probes$START <- ifelse(probes$RSTART < probes$REND, probes$RSTART, probes$REND)
probes$END <- ifelse(probes$RSTART < probes$REND, probes$REND, probes$RSTART)

which(probes$END < probes$START) # 0

probes <- probes[ , c(6, 11, 12, 9, 2)]

write.table(probes, "CornellHempDesign1v3_IntList.txt", quote = F,
            row.names = F, col.names = F, sep = "\t")

final <- merge(final, probes, by.x = "ProbeID", by.y = "PRBID", all = T)

sum(final$Replication) # 54201

lapply(unique(final$CHR), function(x) {

  c(x, sum(final[which(final$CHR == x), ]$REP))

  })
