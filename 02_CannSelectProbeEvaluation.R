# CannSelect probe evaluation

library(stringr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(scales)
library(qqman)
library(seqinr)

GenomeBlastResults <- function(path, version) {

  prbct <- ifelse(version == "V2", 19452, 18040)

  output <- NULL

  resfiles <- list.files(path = path, pattern = ".txt")

  filect <- length(resfiles)
  ct <- 0
  thres <- c(1, 5, 25, 50, 75, 90, 95, 100)
  nxtthrs <- 1

  for (bres in resfiles) {

    res <- read.table(paste(path, bres, sep = "/"))
    res$JScore <- (res$V4 - (res$V5 + res$V6))

    rawct <- nrow(res)
    assembly <- word(bres, 1, sep = "[.]")

    res <- droplevels(res[which(res$JScore >= 60), ])

    aligncounts <- as.data.frame(table(res$V1))

    aligned <- nrow(aligncounts)
    percalign <- round(((aligned/prbct)*100), 1)

    aligncounts$SPEC <- ifelse(aligncounts$Freq < 3, "YES", "NO")

    specific <- sum(aligncounts$SPEC == "YES")
    percspec <- round(((specific/prbct)*100), 1)

    outrow <- c(assembly, rawct, aligned, percalign, specific, percspec)

    output <- rbind(output, outrow)

    ct <- ct + 1

    percdone <- (ct/filect)*100
    if (percdone >= thres[nxtthrs]) {
      message(sprintf("Processed file %d of %d (%.1f%% Complete)",
                      ct, filect, percdone))
      nxtthrs <- nxtthrs + 1
    }
  }

  output <- as.data.frame(output)
  rownames(output) <- NULL
  colnames(output) <- c("Assembly", "Raw",
                        "Aligned", "Perc.Aligned",
                        "Specific", "Perc.Specific")
  colnames(output)[2:6] <- paste(colnames(output)[2:6], version, sep = "_")

  return(output)
}

GetDist2NextProbe <- function(design) {

  output <- NULL

  yprobes <- droplevels(design[which(is.na(design$REFCHR)), ])
  design <- droplevels(design[which(!(design$ProbeID %in% yprobes$ProbeID)), ])

  for (chrom in unique(design$REFCHR)) {

    chrout <- droplevels(design[which(design$REFCHR == chrom), ])
    chrout <- chrout[order(chrout$START) , ]

    chrout$Dist2Next <- NA
    for (prb in seq_len(nrow(chrout))) {
      if (prb == nrow(chrout)) {
        chrout$Dist2Next[prb] <- 0
      } else {
        chrout$Dist2Next[prb] <- chrout$START[prb + 1] - chrout$STOP[prb]
      }
    }
    chrout$Dist2Next <- ifelse(chrout$Dist2Next <= 0, NA, chrout$Dist2Next)

    output <- rbind(output, chrout)
  }

  yprobes$Dist2Next <- NA
  output <- rbind(output, yprobes)

  return(output)
}

GetBinnedProbeCounts <- function(chrtable, binsize, design) {

  design <- droplevels(design[which(!(is.na(design$REFCHR))), ])

  bins <- NULL
  for (chrom in chrtable$CHR) {
    binct <- ceiling(chrtable$SIZE[which(chrtable$CHR == chrom)] / binsize)
    chrout <- as.data.frame(list("CHR" = rep(chrom, binct),
                                 "START" = (seq(1, binct*binsize, by = binsize)-1),
                                 "END" = (seq((1e6-1), binct*binsize, by = binsize))))
    bins <- rbind(bins, chrout)
  }

  bins$PRBCT <- unlist(lapply(seq_len(nrow(bins)), function(x) {
    sum(design$REFCHR == bins$CHR[x] & design$MidPoint >= bins$START[x] & design$MidPoint <= bins$END[x])
  }))

  return(bins)
}

# Cross-assembly specificity
###### 1-2. Getting Summary Results Table ######

path <- "PangenomeBlast/blast_results_probesetv2"
v2 <- GenomeBlastResults(path, "V2")

path <- "PangenomeBlast/blast_results_probesetv3"
v3 <- GenomeBlastResults(path, "V3")

res <- merge(v2, v3, by = "Assembly", all = T)

divgroups <- read.csv("PangenomeBlast/CannSelect_PangenomeDivGroups.csv")

meta <- merge(divgroups, res, by = "Assembly", all = T)

write.csv(meta, "PangenomeBlast/CannSelect_AssemblySpecificityTable.csv",
          row.names = F, quote = F)

# reading it back in
meta <- read.csv("PangenomeBlast/CannSelect_AssemblySpecificityTable.csv")
meta$Group <- as.factor(meta$Group)

meta$Group <- ifelse(meta$Group == "OCBD F1", "High-Cannabinoid", meta$Group)

euros <- c("FCS1b.softmasked", "H3S1a.softmasked",
           "H3S1b.softmasked", "H3S7b.softmasked")
meta$Group <- ifelse(meta$FullAssemblyName %in% euros, "European Hemp", meta$Group)

# getting group means and standard error

divres_mean <- aggregate.data.frame(x = meta[ , c(7, 9, 12, 14)], by = list("GROUP" =  meta$Group), FUN = mean)

divres_se <- aggregate.data.frame(x = meta[ , c(7, 9, 12, 14)], by = list("GROUP" =  meta$Group), FUN = function(x) sd(x) / sqrt(length(x)))

divres_min <- aggregate.data.frame(x = meta[ , c(7, 9, 12, 14)], by = list("GROUP" =  meta$Group), FUN = min)
divres_max <- aggregate.data.frame(x = meta[ , c(7, 9, 12, 14)], by = list("GROUP" =  meta$Group), FUN = max)

colnames(divres_se)[-1] <- paste(colnames(divres_se)[-1], "_SE", sep = "")
colnames(divres_min)[-1] <- paste(colnames(divres_min)[-1], "_MIN", sep = "")
colnames(divres_max)[-1] <- paste(colnames(divres_max)[-1], "_MAX", sep = "")

divres <- merge(divres_mean,
                merge(divres_se,
                      merge(divres_min, divres_max, by = "GROUP"),
                      by = "GROUP"),
                by = "GROUP")

write.csv(divres, "PangenomeBlast/CannSelect_GroupSpecificityTable_4Groups.csv",
          row.names = F, quote = F)

# CannSelect v0 capture metrics and figures
###### 3-1. Sample and Target Statistics ######

### SAMPLE METRICS ###

path <- "HSMET"

samplist <- list.files(path, pattern = "HSMET")
samplist <- samplist[which(grepl("CU.HSMET", samplist))]

sampmet <- read.table(paste(path, samplist[1], sep = "/"),
                      sep = "\t", skip = 6, nrows = 2) [ , -c(1, 67:69)]
colnames(sampmet) <- sampmet[ 1, ]
sampmet <- sampmet[-1, ]
sampmet$SAMPID <- word(samplist[1], 1, sep = "[.]")
sampmet <- sampmet[ , c(66, 1:65)]

samphist <- read.table(paste(path, samplist[1], sep = "/"),
                       sep = "\t", skip = 10, header = T)[ , 1:2]
samphist$SAMPID <- word(samplist[1], 1, sep = "[.]")
colnames(samphist) <- c("COV", "HQCT", "SAMPID")

for (samp in samplist[-1]) {

  met <- read.table(paste(path, samp, sep = "/"),
                    sep = "\t", skip = 6, nrows = 2) [ , -c(1, 67:69)]
  colnames(met) <- met[ 1, ]
  met <- met[-1, ]
  met$SAMPID <- word(samp, 1, sep = "[.]")
  met <- met[ , c(66, 1:65)]

  hst <- read.table(paste(path, samp, sep = "/"),
                    sep = "\t", skip = 10, header = T)[ , 1:2]
  hst$SAMPID <- word(samp, 1, sep = "[.]")
  colnames(hst) <- c("COV", "HQCT", "SAMPID")

  sampmet <- rbind(sampmet, met)
  samphist <- rbind(samphist, hst)

}

hist(as.numeric(sampmet$PCT_SELECTED_BASES))
write.table(samphist, "SampHSMet_CovHist.txt", quote = F, row.names = F)
write.table(sampmet, "SampHSMet_Stats.txt", quote = F, row.names = F)

### INTERVAL METRICS ###

intlist <- list.files(path, pattern = "PERTARG")
intlist <- intlist[which(grepl("CU.PERTARG", intlist))]

prblist <- read.table("CornellHempDesign1v2_Probes.txt", header = T, sep = "\t")[ , c(1:4)]

intmet <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:6)]

intmet$PRBCT <- unlist(lapply(1:nrow(intmet), function(x) {
  length(str_split(intmet$name[x], "[|]")[[1]])
}))

intmet$REPCT <- unlist(lapply(1:nrow(intmet), function(x) {
  if (intmet$PRBCT[x] == 1) {
    prblist$Replication[which(prblist$ProbeID == intmet$name[x])]
  } else {
    count <- 0
    for (pro in str_split(intmet$name[x], "[|]")[[1]]) {
      count <- count + prblist$Replication[which(prblist$ProbeID == pro)]
    }
    count
  }
}))

intmet$name <- NULL

avgcov <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:4, 6, 7)]
mincov <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:4, 6, 11)]
maxcov <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:4, 6, 12)]
p0x <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:4, 6, 13)]
rct <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:4, 6, 14)]
first <- c("REFCHR", "START", "END", "LENGTH", "PGC", word(intlist[1], 1, sep = "[.]"))
colnames(intmet)[1:5] <- first[1:5]
colnames(avgcov) <- first
colnames(mincov) <- first
colnames(maxcov) <- first
colnames(p0x) <- first
colnames(rct) <- first

for (int in intlist[-1]) {
  met <- read.table(paste(path, int, sep = "/"), header = T)[ , c(1:4, 6, 7, 11:14)]
  colnames(met)[1:5] <- first[1:5]
  avgcov <- merge(avgcov, met[ , c(1:5, 6)], by = first[1:5], all = T)
  mincov <- merge(mincov, met[ , c(1:5, 7)], by = first[1:5], all = T)
  maxcov <- merge(maxcov, met[ , c(1:5, 8)], by = first[1:5], all = T)
  p0x <- merge(p0x, met[ , c(1:5, 9)], by = first[1:5], all = T)
  rct <- merge(rct, met[ , c(1:5, 10)], by = first[1:5], all = T)
  colnames(avgcov)[ncol(avgcov)] <- word(int,  1, sep = "[.]")
  colnames(mincov)[ncol(mincov)] <- word(int,  1, sep = "[.]")
  colnames(maxcov)[ncol(maxcov)] <- word(int,  1, sep = "[.]")
  colnames(p0x)[ncol(p0x)] <- word(int,  1, sep = "[.]")
  colnames(rct)[ncol(rct)] <- word(int,  1, sep = "[.]")
}

avgcov <- merge(intmet, avgcov, by = first[1:5], all = T)
mincov <- merge(intmet, mincov, by = first[1:5], all = T)
maxcov <- merge(intmet, maxcov, by = first[1:5], all = T)
p0x <- merge(intmet, p0x, by = first[1:5], all = T)
rct <- merge(intmet, rct, by = first[1:5], all = T)

write.table(avgcov, "IntHSMet_MeanCov.txt", quote = F, sep = "\t", row.names = F)
write.table(mincov, "IntHSMet_MinCov.txt", quote = F, sep = "\t", row.names = F)
write.table(maxcov, "IntHSMet_MaxCov.txt", quote = F, sep = "\t", row.names = F)
write.table(p0x, "IntHSMet_PZeroCov.txt", quote = F, sep = "\t", row.names = F)
write.table(rct, "IntHSMet_ReadCount.txt", quote = F, sep = "\t", row.names = F)

###### 4-1. CannSelect.v2 Run Figure ######

# 6in x 6in @ 300 dpi = 1800x1800 px

### coverage by chromosome ###

meancov <- read.table("IntHSMet_MeanCov.txt", header = T)
meancov$SAMPMEAN <- as.numeric(unlist(lapply(seq_len(nrow(meancov)), function(x) {
  mean(as.numeric(meancov[x, c(8:ncol(meancov))]))
})))

meancov$CHR <- ifelse(meancov$REFCHR == "NC_044371.1", "CHR01",
                      ifelse(meancov$REFCHR == "NC_044375.1", "CHR02",
                             ifelse(meancov$REFCHR == "NC_044372.1", "CHR03",
                                    ifelse(meancov$REFCHR == "NC_044373.1", "CHR04",
                                           ifelse(meancov$REFCHR == "NC_044374.1", "CHR05",
                      ifelse(meancov$REFCHR == "NC_044377.1", "CHR06",
                             ifelse(meancov$REFCHR == "NC_044378.1", "CHR07",
                                    ifelse(meancov$REFCHR == "NC_044379.1", "CHR08",
                                           ifelse(meancov$REFCHR == "NC_044376.1", "CHR09",
                                                  ifelse(meancov$REFCHR == "NC_044370.1", "CHRX", "SCAFFOLD"))))))))))

meancov$CHR <- factor(meancov$CHR, levels = c("CHR01", "CHR02", "CHR03", "CHR04", "CHR05",
                                              "CHR06", "CHR07", "CHR08", "CHR09", "CHRX"))

covsummary <- meancov %>%
  group_by(CHR) %>%
  summarise(CHRMEAN = mean(SAMPMEAN),
            SE = sd(SAMPMEAN) / sqrt(n()))

covbychr <- ggplot(data = covsummary, aes(x = CHR, y = CHRMEAN)) +
#ggplot(data = covsummary, aes(x = CHR, y = CHRMEAN)) +
  geom_col(fill = "gray75") +
  geom_errorbar(aes(ymin = CHRMEAN - SE, ymax = CHRMEAN + SE),
                col = "firebrick2", width = 0.5, linewidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.25, 42)) +
  xlab("Chromosome (CBDRx/cs10)") +
  ylab("Mean Target Fold Coverage") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(fill = NA),
        axis.title = element_text(size = 8, face = "bold", color = "black"),
        axis.text.x = element_text(size = 7, face = "bold", color = "black",
                                   angle = 45, vjust = .6),
        axis.text.y = element_text(size = 7, face = "bold", color = "black"))

ggsave(filename = "Figures/TargetCoverageByChr.pdf",
       plot = covbychr,
       device = cairo_pdf,
       width = 3, height = 3,
       units = "in",
       dpi = 350)

### Sample Metrics ###

samp <- read.table("SampHSMet_Stats.txt", header = T)

# this one is dual axis...
count_cols <- c("PF_READS", "PF_UNIQUE_READS", "PF_UQ_READS_ALIGNED")
pct_cols   <- c("PCT_SELECTED_BASES", "ZERO_CVG_TARGETS_PCT")

df <- samp

# --- 2) Compute a scaling factor to align % with counts on the left axis ---
# Convert percent columns to 0–100 if they are in 0–1; otherwise leave as-is
pct_vals_raw <- unlist(df[pct_cols])
pct_vals_0100 <- if (max(pct_vals_raw, na.rm = TRUE) <= 1) pct_vals_raw * 100 else pct_vals_raw

max_count <- 4.7e6
max_pct   <- 100
s <- if (is.finite(max_count) && is.finite(max_pct) && max_pct > 0) max_count / max_pct else 1

# --- 3) Reshape long and create plotting values with scaling applied to % ---
df_long <-
  bind_rows(
    df[count_cols] |>
      pivot_longer(everything(), names_to = "var", values_to = "value") |>
      mutate(type = "count", value_plot = value),
    # percent side (convert to 0–100, then scale by s so it fits the count axis)
    df[pct_cols] |>
      pivot_longer(everything(), names_to = "var", values_to = "value") |>
      mutate(type = "pct",
             value_0100 = if (max(value, na.rm = TRUE) <= 1) value * 100 else value,
             value_plot = value_0100 * s)
  ) |>
  select(var, type, value, value_plot)

# --- 4) Summary stats (mean ± SE) on the plotting scale ---
summ_df <- df_long |>
  group_by(var, type) |>
  summarise(
    n = sum(!is.na(value_plot)),
    mean_plot = mean(value_plot, na.rm = TRUE),
    se_plot   = sd(value_plot, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

summ_df$var <- factor(summ_df$var, levels = c("PF_READS", "PF_UNIQUE_READS", "PF_UQ_READS_ALIGNED",
                                              "PCT_SELECTED_BASES", "ZERO_CVG_TARGETS_PCT"))
df_long$var <- factor(df_long$var, levels = c("PF_READS", "PF_UNIQUE_READS", "PF_UQ_READS_ALIGNED",
                                              "PCT_SELECTED_BASES", "ZERO_CVG_TARGETS_PCT"))

sampdat <- ggplot() +
#ggplot() +
  geom_col(data = summ_df, aes(var, mean_plot), fill = "gray75") +
  geom_jitter(data = df_long, aes(var, value_plot), shape = 19,
              width = 0.3, height = 0, alpha = 0.2, size = 1) +
  geom_errorbar(data = summ_df, aes(var, ymin = mean_plot - se_plot, ymax = mean_plot + se_plot),
                width = 0.5, color = "firebrick2", linewidth = 0.4) +
  geom_vline(aes(xintercept = 3.5), linewidth = 0.1) +
  scale_y_continuous(name = "Read Count (Millions)", sec.axis = sec_axis(~ . / s, name = "Percent"),
                     expand = c(0.01, 0.01), limits = c(-500, 4.7e6),
                     labels = seq.int(0, 4, by = 1),
                     breaks = seq.int(0, 4e6, by = 1e6)) +
  scale_x_discrete(labels = c("Total\nReads", "Unique\nReads", "Aligned\nReads",
                              "Percent\nSelected", "Percent\n Zero Cov."),
                   expand = c(0.11, 0.11)) +
  labs(x = NULL) +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(fill = NA),
        axis.title = element_text(size = 8, face = "bold", color = "black"),
        axis.text.x = element_text(size = 7, face = "bold", color = "black"),
        axis.text.y = element_text(size = 7, face = "bold", color = "black"))

ggsave(filename = "Figures/v2SampleMetrics.pdf",
       plot = sampdat,
       device = cairo_pdf,
       width = 3, height = 3,
       units = "in",
       dpi = 350)

# Design statistics and density plots
###### 4-2. CannSelect Design Statistics ######

v2probes <- read.table("CornellHempDesign1v2_Probes.txt", header = T, sep = "\t")
v3probes <- read.table("CornellHempDesign1v3_Probes.txt", header = T, sep = "\t")
v2ints <- read.table("CornellHempDesign1v2_IntList.txt")
v3ints <- read.table("CornellHempDesign1v3_IntList.txt")
colnames(v2ints) <- colnames(v3ints) <- c("REFCHR", "START", "STOP", "STRAND", "ProbeID")

v3ints$ProbeID <- str_replace(v3ints$ProbeID, "_ProbeID", "")

csv2 <- merge(v2probes[ , c(1:2, 4)], v2ints[ , -c(4)], by = "ProbeID", all = T)
csv3 <- merge(v3probes[ , c(1:2, 4)], v3ints[ , -c(4)], by = "ProbeID", all = T)

rm(v2probes, v2ints, v3probes, v3ints)

# setting marker v gene based
csv2$ProbeType <- ifelse(grepl("NC", csv2$TargetID), "Marker-Based", "Gene-Based")
csv3$ProbeType <- ifelse(grepl("NC", csv3$TargetID), "Marker-Based", "Gene-Based")

# determining probe mid point for plotting
csv2$MidPoint <- (csv2$START + csv2$STOP) / 2
csv3$MidPoint <- (csv3$START + csv3$STOP) / 2

# finding distance to next probe (sets negatives to missing)
csv2 <- GetDist2NextProbe(csv2)
csv3 <- GetDist2NextProbe(csv3)

# finding probes per Mb  #

chrtable <- as.data.frame(list("CHR" = c("NC_044371.1", "NC_044375.1", "NC_044372.1",
                                         "NC_044373.1", "NC_044374.1", "NC_044377.1",
                                         "NC_044378.1", "NC_044379.1", "NC_044376.1",
                                         "NC_044370.1"),
                               "SIZE" = c(101209240, 96346938, 94670641,
                                          91913879, 88181582, 79335105,
                                          71238074, 64622176, 61561104,
                                          104987320)))

csv2_bins <- GetBinnedProbeCounts(chrtable, 1e6, csv2)
csv3_bins <- GetBinnedProbeCounts(chrtable, 1e6, csv3)

### getting stats:v2 ###

# probeset total = 19452

# probe type
table(csv2$ProbeType)
# Gene-Based Marker-Based

# boosted total
sum(csv2$Replication) # 54201

# probes per chrom
table(csv2$REFCHR)
# NC_044370.1 NC_044371.1 NC_044372.1 NC_044373.1 NC_044374.1
# NC_044375.1 NC_044376.1 NC_044377.1 NC_044378.1 NC_044379.1

# probes by type by chrom
table(csv2$REFCHR, csv2$ProbeType)
#             Gene-Based Marker-Based
# NC_044370.1        328         2191
# NC_044371.1        281         2057
# NC_044372.1        112         1606
# NC_044373.1        183         1848
# NC_044374.1        105         1799
# NC_044375.1        301         1869
# NC_044376.1        335         1358
# NC_044377.1        119         1568
# NC_044378.1        175         1296
# NC_044379.1        369         1524

# boosted total per chromosome
boost <- unlist(lapply(unique(csv2$REFCHR), function(x) {
  if (!is.na(x)) {
    sum(csv2$Replication[which(csv2$REFCHR == x)])
  } else {
    sum(csv2$Replication[which(is.na(csv2$REFCHR))])
  }
}))
names(boost) <- unique(csv2$REFCHR)
boost
# NC_044370.1 NC_044371.1 NC_044372.1 NC_044373.1 NC_044374.1
# NC_044375.1 NC_044376.1 NC_044377.1 NC_044378.1 NC_044379.1
# <NA>

# avg probe count per mb by chr
prbct <- unlist(lapply(unique(csv2_bins$CHR), function(x) {
    mean(csv2_bins$PRBCT[which(csv2_bins$CHR == x)])
}))
names(prbct) <- unique(csv2_bins$CHR)
prbct
# NC_044371.1 NC_044375.1 NC_044372.1 NC_044373.1 NC_044374.1
# NC_044377.1 NC_044378.1 NC_044379.1 NC_044376.1 NC_044370.1

# maximum gap size per chr
maxgap <- unlist(lapply(unique(csv2$REFCHR), function(x) {
  if (!is.na(x)) {
    max(csv2$Dist2Next[which(csv2$REFCHR == x)], na.rm = T)
  } else { NA }
}))
names(maxgap) <- unique(csv2$REFCHR)
maxgap/1000000 # div so nums are in Mb
# NC_044370.1 NC_044371.1 NC_044372.1 NC_044373.1 NC_044374.1
# NC_044375.1 NC_044376.1 NC_044377.1 NC_044378.1 NC_044379.1
# <NA>
#   NA

### getting stats:v3 ###

# probeset total = 18040

# probe type
table(csv3$ProbeType)
# Gene-Based Marker-Based

# boosted total
sum(csv3$Replication) # 54201

# probes per chrom
table(csv3$REFCHR, useNA = "ifany")
# NC_044370.1 NC_044371.1 NC_044372.1 NC_044373.1 NC_044374.1
# NC_044375.1 NC_044376.1 NC_044377.1 NC_044378.1 NC_044379.1
# <NA>

# probes by type by chrom
table(csv3$REFCHR, csv3$ProbeType)
#             Gene-Based Marker-Based
# NC_044370.1        417         1843
# NC_044371.1        392         1889
# NC_044372.1        168         1468
# NC_044373.1        225         1644
# NC_044374.1        138         1495
# NC_044375.1        322         1639
# NC_044376.1        362         1260
# NC_044377.1        128         1372
# NC_044378.1        187         1104
# NC_044379.1        539         1420

# boosted total per chromosome
boost <- unlist(lapply(unique(csv3$REFCHR), function(x) {
  if (!is.na(x)) {
    sum(csv3$Replication[which(csv3$REFCHR == x)])
  } else {
    sum(csv3$Replication[which(is.na(csv3$REFCHR))])
  }
}))
names(boost) <- unique(csv3$REFCHR)
boost
# NC_044378.1 NC_044370.1 NC_044371.1 NC_044372.1 NC_044373.1
# NC_044374.1 NC_044375.1 NC_044376.1 NC_044377.1 NC_044379.1
# <NA>

# avg probe count per mb by chr
prbct <- unlist(lapply(unique(csv3_bins$CHR), function(x) {
  mean(csv3_bins$PRBCT[which(csv3_bins$CHR == x)])
}))
names(prbct) <- unique(csv3_bins$CHR)
prbct
# NC_044371.1 NC_044375.1 NC_044372.1 NC_044373.1 NC_044374.1
# NC_044377.1 NC_044378.1 NC_044379.1 NC_044376.1 NC_044370.1

# maximum gap size per chr
maxgap <- unlist(lapply(unique(csv3$REFCHR), function(x) {
  if (!is.na(x)) {
    max(csv3$Dist2Next[which(csv3$REFCHR == x)], na.rm = T)
  } else { NA }
}))
names(maxgap) <- unique(csv3$REFCHR)
maxgap/1000000 # div so nums are in Mb
# NC_044378.1 NC_044370.1 NC_044371.1 NC_044372.1 NC_044373.1
# NC_044374.1 NC_044375.1 NC_044376.1 NC_044377.1 NC_044379.1
# <NA>
#   NA
###### 4-3. CannSelect.v3 Density Plot ######

v3probes <- read.table("CornellHempDesign1v3_Probes.txt", header = T, sep = "\t")
v3ints <- read.table("CornellHempDesign1v3_IntList.txt")
colnames(v3ints) <- c("REFCHR", "START", "STOP", "STRAND", "ProbeID")

v3ints$ProbeID <- str_replace(v3ints$ProbeID, "_ProbeID", "")
csv3 <- merge(v3probes[ , c(1:2, 4)], v3ints[ , -c(4)], by = "ProbeID", all = T)

#rm(v3probes, v3ints)

# setting marker v gene based
csv3$ProbeType <- ifelse(grepl("NC", csv3$TargetID), "Marker-Based", "Gene-Based")

# determining probe mid point for plotting
csv3$MidPoint <- (csv3$START + csv3$STOP) / 2

# finding probes per Mb  #

chrtable <- as.data.frame(list("CHR" = c("NC_044371.1", "NC_044375.1", "NC_044372.1",
                                         "NC_044373.1", "NC_044374.1", "NC_044377.1",
                                         "NC_044378.1", "NC_044379.1", "NC_044376.1",
                                         "NC_044370.1"),
                               "SIZE" = c(101209240, 96346938, 94670641,
                                          91913879, 88181582, 79335105,
                                          71238074, 64622176, 61561104,
                                          104987320)))

csv3_bins <- GetBinnedProbeCounts(chrtable, 1e6, csv3)

# setting chr y positions
chrtable$ymin_den <- seq.int(2, 51, by = 5)
chrtable$ymax_den <- seq.int(5, 51, by = 5)
chrtable$ymin_prb <- seq.int(5, 51, by = 5)
chrtable$ymax_prb <- seq.int(6, 51, by = 5)

csv3 <- merge(csv3[ , 4:8], chrtable, by.x = "REFCHR", by.y = "CHR")
csv3_bins <- merge(csv3_bins, chrtable, by = "CHR")

# coverting to Mb
csv3$plotMB <- round(csv3$MidPoint/1e6, 1)
csv3_bins$STARTMB <- round(csv3_bins$START/1e6, 1)
csv3_bins$ENDMB <- round(csv3_bins$END/1e6, 1)

# setting bin colors
csv3_bins$plotct <- ifelse(csv3_bins$PRBCT >= 50, 50, csv3_bins$PRBCT)
csv3_bins$prbct_norm <- rescale(csv3_bins$plotct, to = c(0, 1))
ramp <- colorRamp(c("white", "gray50", "firebrick3"))
csv3_bins$bincol <- rgb(ramp(csv3_bins$prbct_norm) / 255)

# making the plot
denplot <- ggplot( ) +
#ggplot( ) +
  scale_x_continuous(expand = c(0,0), limits = c(-0.5, 106),
                     breaks = seq.int(0, 105, by = 20),
                     labels = seq.int(0, 105, by = 20)) +
  scale_y_reverse(expand = c(0,0), limits = c(52, 1),
                  breaks = seq.int(4, 52, by = 5),
                  labels = c("CHR01\nNC_044371.1", "CHR02\nNC_044375.1", "CHR03\nNC_044372.1",
                             "CHR04\nNC_044373.1", "CHR05\nNC_044374.1", "CHR06\nNC_044377.1",
                             "CHR07\nNC_044378.1", "CHR08\nNC_044379.1", "CHR09\nNC_044376.1",
                             "CHRX\nNC_044370.1")) +
  geom_rect(data = csv3_bins, aes(xmin = STARTMB, xmax = ENDMB,
                                  ymin = ymin_den, ymax = ymax_den),
            fill = csv3_bins$bincol) +
  geom_segment(data = csv3, aes(x = plotMB, xend = plotMB,
                                y = ymin_prb, yend = ymax_prb),
               linewidth = 0.05) +
  geom_rect(data = chrtable, aes(xmin = rep(0, 10), xmax = ceiling(SIZE/1e6),
                                 ymin = ymin_den, ymax = ymax_prb),
            fill = NA, col = "black", linewidth = 0.25) +
  geom_segment(data = chrtable, aes(x = rep(0, 10), xend =  ceiling(SIZE/1e6),
                                    y = ymax_den, yend = ymax_den),
               linewidth = 0.25) +
  xlab("CBDRx(cs10) Physical Position (Mb)") +
  ylab("CBDRx(cs10) Chromosome") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(fill = NA),
        axis.title = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(size = 7, face = "bold", color = "black"),
        axis.text.y = element_text(size = 7, face = "bold", color = "black"))

ggsave(filename = "Figures/CannSelectv3DesignDensity.pdf",
       plot = denplot,
       device = cairo_pdf,
       width = 5.25, height = 4,
       units = "in",
       dpi = 350)

# making the gradient legend
binlegend <- droplevels(csv3_bins[which(!duplicated(csv3_bins$bincol)), ])
binlegend <- binlegend[order(binlegend$prbct_norm), ]

binlegend$ypos <- seq.int(1, 47)

legend_plot <- ggplot(binlegend, aes(x = 1, y = ypos, fill = bincol)) +
#ggplot(binlegend, aes(x = 1, y = ypos, fill = bincol)) +
  geom_tile() +
  scale_fill_identity() +
  scale_y_continuous(breaks = c(1,
                                (max(binlegend$ypos)/5)*1,
                                (max(binlegend$ypos)/5)*2,
                                (max(binlegend$ypos)/5)*3,
                                (max(binlegend$ypos)/5)*4,
                                max(binlegend$ypos)),
                     labels = c(0, 10, 20, 30, 40, "50+"),
                     position = "right",
                     expand = c(0.001, 0.001),
                     name = NULL) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = NULL) +
  ggtitle("Probes\nPer Mb") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 8, face = "bold",
                             color = "black", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7, face = "bold", color = "black"))

ggsave(filename = "Figures/CannSelectv3DesignDensity_Legend.pdf",
       plot = legend_plot,
       device = cairo_pdf,
       width = 0.75, height = 3.5,
       units = "in",
       dpi = 350)

# CannSelect v1.0 capture metrics and figure
###### 5-3. USDAv3 Sample and Target Statistics ######

### SAMPLE METRICS ###

path <- "USDAv3/HSMET_USDA"

samplist <- list.files(path, pattern = "HSMET")
# removing blanks
samplist <- samplist[which(!grepl("BLANK", samplist))]

sampmet <- read.table(paste(path, samplist[1], sep = "/"),
                      sep = "\t", skip = 6, nrows = 2) [ , -c(1, 67:69)]
colnames(sampmet) <- sampmet[ 1, ]
sampmet <- sampmet[-1, ]
sampmet$SAMPID <- word(samplist[1], 1, sep = "[.]")
sampmet <- sampmet[ , c(66, 1:65)]

samphist <- read.table(paste(path, samplist[1], sep = "/"),
                       sep = "\t", skip = 10, header = T)[ , 1:2]
samphist$SAMPID <- word(samplist[1], 1, sep = "[.]")
colnames(samphist) <- c("COV", "HQCT", "SAMPID")

for (samp in samplist[-1]) {

  met <- read.table(paste(path, samp, sep = "/"),
                    sep = "\t", skip = 6, nrows = 2) [ , -c(1, 67:69)]
  colnames(met) <- met[ 1, ]
  met <- met[-1, ]
  met$SAMPID <- word(samp, 1, sep = "[.]")
  met <- met[ , c(66, 1:65)]

  hst <- read.table(paste(path, samp, sep = "/"),
                    sep = "\t", skip = 10, header = T)[ , 1:2]
  hst$SAMPID <- word(samp, 1, sep = "[.]")
  colnames(hst) <- c("COV", "HQCT", "SAMPID")

  sampmet <- rbind(sampmet, met)
  samphist <- rbind(samphist, hst)

}

mean(as.numeric(sampmet$PCT_SELECTED_BASES)) # 88.2
write.table(samphist, "USDAv3/SampHSMet_CovHist.txt", quote = F, row.names = F)
write.table(sampmet, "USDAv3/SampHSMet_Stats.txt", quote = F, row.names = F)

### INTERVAL METRICS ###

intlist <- list.files(path, pattern = "PERTARG")
intlist <- intlist[which(!grepl("BLANK", intlist))]

prblist <- read.table("CornellHempDesign1v3_Probes.txt", header = T, sep = "\t")[ , c(1:4)]

intmet <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:6)]

intmet$PRBCT <- unlist(lapply(1:nrow(intmet), function(x) {
  length(str_split(intmet$name[x], "[|]")[[1]])
}))

intmet$name <- str_replace_all(intmet$name, "ProbeID_", "")

intmet$REPCT <- unlist(lapply(1:nrow(intmet), function(x) {
  if (intmet$PRBCT[x] == 1) {
    prblist$Replication[which(prblist$ProbeID == intmet$name[x])]
  } else {
    count <- 0
    for (pro in str_split(intmet$name[x], "[|]")[[1]]) {
      count <- count + prblist$Replication[which(prblist$ProbeID == pro)]
    }
    count
  }
}))

intmet$name <- NULL

avgcov <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:4, 6, 7)]
mincov <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:4, 6, 11)]
maxcov <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:4, 6, 12)]
p0x <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:4, 6, 13)]
rct <- read.table(paste(path, intlist[1], sep = "/"), header = T)[ , c(1:4, 6, 14)]
first <- c("REFCHR", "START", "END", "LENGTH", "PGC", word(intlist[1], 1, sep = "[.]"))
colnames(intmet)[1:5] <- first[1:5]
colnames(avgcov) <- first
colnames(mincov) <- first
colnames(maxcov) <- first
colnames(p0x) <- first
colnames(rct) <- first

for (int in intlist[-1]) {
  met <- read.table(paste(path, int, sep = "/"), header = T)[ , c(1:4, 6, 7, 11:14)]
  colnames(met)[1:5] <- first[1:5]
  avgcov <- merge(avgcov, met[ , c(1:5, 6)], by = first[1:5], all = T)
  mincov <- merge(mincov, met[ , c(1:5, 7)], by = first[1:5], all = T)
  maxcov <- merge(maxcov, met[ , c(1:5, 8)], by = first[1:5], all = T)
  p0x <- merge(p0x, met[ , c(1:5, 9)], by = first[1:5], all = T)
  rct <- merge(rct, met[ , c(1:5, 10)], by = first[1:5], all = T)
  colnames(avgcov)[ncol(avgcov)] <- word(int,  1, sep = "[.]")
  colnames(mincov)[ncol(mincov)] <- word(int,  1, sep = "[.]")
  colnames(maxcov)[ncol(maxcov)] <- word(int,  1, sep = "[.]")
  colnames(p0x)[ncol(p0x)] <- word(int,  1, sep = "[.]")
  colnames(rct)[ncol(rct)] <- word(int,  1, sep = "[.]")
}

avgcov <- merge(intmet, avgcov, by = first[1:5], all = T)
mincov <- merge(intmet, mincov, by = first[1:5], all = T)
maxcov <- merge(intmet, maxcov, by = first[1:5], all = T)
p0x <- merge(intmet, p0x, by = first[1:5], all = T)
rct <- merge(intmet, rct, by = first[1:5], all = T)

write.table(avgcov, "USDAv3/IntHSMet_MeanCov.txt", quote = F, sep = "\t", row.names = F)
write.table(mincov, "USDAv3/IntHSMet_MinCov.txt", quote = F, sep = "\t", row.names = F)
write.table(maxcov, "USDAv3/IntHSMet_MaxCov.txt", quote = F, sep = "\t", row.names = F)
write.table(p0x, "USDAv3/IntHSMet_PZeroCov.txt", quote = F, sep = "\t", row.names = F)
write.table(rct, "USDAv3/IntHSMet_ReadCount.txt", quote = F, sep = "\t", row.names = F)

###### 5-4. USDAv3 Figure ######

sampstats <- read.table("USDAv3/SampHSMet_Stats.txt", header = T)

# renaming taxa
sampstats$SAMPID <- word(sampstats$SAMPID, 3, -4, sep = "_")
sampstats$SAMPID <- ifelse(sampstats$SAMPID == "5", "No.5", sampstats$SAMPID)

# removing unreplicated taxa
duplist <- unique(as.character(sampstats$SAMPID[which(duplicated(as.character(sampstats$SAMPID)))]))
sampstats <- droplevels(sampstats[which(sampstats$SAMPID %in% duplist), ]) # 90

sampstats <- droplevels(sampstats[which(!sampstats$SAMPID %in% c("1v1_THCVXSyrian", "BoAx_akaBaOX",
                                                                 "Hohenthurmer_G_R")), ])

selectsum <- sampstats %>%
  group_by(SAMPID) %>%
  summarise(SELMEAN = mean(PCT_SELECTED_BASES),
            SE = sd(PCT_SELECTED_BASES) / sqrt(n()))

pselbytaxa <- ggplot(data = selectsum, aes(x = SAMPID, y = (SELMEAN * 100))) +
  #ggplot(data = selectsum, aes(x = SAMPID, y = (SELMEAN * 100))) +
  geom_col(fill = "gray75") +
  geom_jitter(data = sampstats, aes(x = SAMPID, y = (PCT_SELECTED_BASES*100)),
              shape = 19, width = 0.1, height = 0, alpha = 0.5, size = 1) +
  geom_errorbar(aes(ymin = (SELMEAN - SE)*100, ymax = (SELMEAN + SE)*100),
                col = "firebrick2", width = 0.5, linewidth = 0.5) +
  geom_hline(yintercept = 88.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 101),
                     breaks = c(0, 25, 50, 75, 88.5, 100)) +
  scale_x_discrete(labels = c("Acer II", "Fasamo", "Fibridia", "Goga",
                              "Golden\nRedwood", "Henola", "Iowa Feral 22", "Neo",
                              "No.5", "Picolo", "SS Beta", "Santhica 27",
                              "Tiborszallasi", "YunMa", "Zaki")) +
  xlab("Sample ID") +
  ylab("Percent Selected (%)") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(fill = NA),
        axis.title = element_text(size = 8, face = "bold", color = "black"),
        axis.text.x = element_text(size = 7, face = "bold", color = "black",
                                   angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7, face = "bold", color = "black"))

ggsave(filename = "Figures/USDAv3_PSelectbyTaxa.pdf",
       plot = pselbytaxa,
       device = cairo_pdf,
       width = 6, height = 3,
       units = "in",
       dpi = 350)

# Export final probe sequences to FASTA
### probe sequences to fasta
ver2 <- read.delim("CornellHempDesign1v2_Probes.txt")
ver3 <- read.delim("CornellHempDesign1v3_Probes.txt")

ver2seq <- lapply(ver2$Sequence, function(x) {s2c(x)})
ver3seq <- lapply(ver3$Sequence, function(x) {s2c(x)})

names(ver2seq) <- ver2$ProbeID
names(ver3seq) <- ver3$ProbeID

write.fasta(ver2seq, names(ver2seq),
            file.out = "CannSelectProbeset_v2seq.fasta")
write.fasta(ver3seq, names(ver3seq),
            file.out = "CannSelectProbeset_v3seq.fasta")
