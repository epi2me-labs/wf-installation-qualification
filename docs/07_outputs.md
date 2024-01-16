Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | ./wf-alignment-report.html | Report for all samples | aggregated |
| Combined references | ./combined-refs.fasta | FASTA file containing all input references. | aggregated |
| Combined references index | ./combined-refs.fasta.fai | Index file for combined references FASTA. | aggregated |
| Per-read alignment stats | ./{{ alias }}.readstats.tsv | Bamstats per-read output TSV file. | per-sample |
| Per-reference alignment stats | ./{{ alias }}.flagstat.tsv | Bamstats flagstat output TSV file. | per-sample |
| Alignments BAM file | ./{{ alias }}.sorted.aligned.bam | BAM file with alignments of filtered input reads against the combined references. | per-sample |
| Alignments index file | ./{{ alias }}.sorted.aligned.bam.bai | Index for alignments BAM file. | per-sample |
