This workflow enables users to carry out a series of tests on each sample, to compare statistics to predefined criteria and thresholds.
These include mapping statistics based on alignment to a specified reference sequence, for both test and negative control samples.
A set of default criteria is provided in data/tests_config.json.
An alignment report is also produced to illustrate and summarise the results.

In brief, the workflow performs the following:
* Combine all reference files in the directory passed to `--references`.
* Align input reads (passed as FASTQ or unaligned BAM files) against the reference (note that BAM files with aligned reads can be used as well; these will skip the alignment step and only statistics and the report will be produced).
* Create alignment statistics.
* Compare alignment statistics to provided criteria and marks samples as failed or passed.
* Create an HTML report to illustrate the results.
