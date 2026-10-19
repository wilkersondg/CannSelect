#!/usr/bin/env bash
# BLAST evaluation of the completed CannSelect probe designs.
# 
set -euo pipefail

cd PangenomeBlast
mkdir -p RefResults

makeblastdb -in RefGenomes/CS10_genomic.fasta -dbtype nucl
makeblastdb -in RefGenomes/PinkPepper_genomic.fasta -dbtype nucl

for version in v2 v3; do
  query="CannSelectProbeset_${version}seq.fasta"
  blastn -db RefGenomes/CS10_genomic.fasta -query "$query" -outfmt 6 -max_target_seqs 5 -evalue 1e-10 -out "RefResults/CS10_${version}.txt"
  blastn -db RefGenomes/PinkPepper_genomic.fasta -query "$query" -outfmt 6 -max_target_seqs 5 -evalue 1e-10 -out "RefResults/PinkPepper_${version}.txt"
done

# Additional pangenome BLAST outputs are summarized by
# GenomeBlastResults() in 02_CannSelectProbeEvaluation.R.
