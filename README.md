# CannSelect probe design and evaluation

Source code accompanying the CannSelect manuscript. This repository contains the custom probe-design and probe-evaluation code only; standard sequencing and variant-calling pipelines are outside its scope.

- `01_CannSelectProbeDesign.R` — gene- and marker-based probe generation, selection, design integration, and final export.
- `02_CannSelectProbeEvaluation.R` — design statistics, genome-wide density, capture metrics, and cross-assembly BLAST summaries.
- `03_CannSelectProbeAlignment.sh` — BLAST searches of the completed designs against the CS10 and Pink Pepper reference genomes.

## Inputs and execution

Scripts use relative file paths and expect the original input tables, FASTA sequences, and precomputed intermediate results in the working directory. These inputs are referenced by their original file names in the scripts. The scripts are provided as analytical source code rather than as a self-contained, one-command pipeline.

The probe-design script also consumes independently produced BLAST results; the alignment shell script documents the final-design evaluation searches, not every upstream design search.

## Dependencies

R packages: `stringr`, `ggplot2`, `reshape`, `seqinr`, `scales`, `dplyr`, `tidyverse`, and `qqman` (where loaded). BLAST+ is required for the alignment shell script.

The code has been consolidated for deposit but has not been rerun against the original inputs.
