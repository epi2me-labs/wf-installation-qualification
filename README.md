# Installation qualification workflow

Align Nanopore reads and visualize mapping statistics.



## Introduction

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




## Compute requirements

Recommended requirements:

+ CPUs = 12
+ Memory = 32GB

Minimum requirements:

+ CPUs = 4
+ Memory = 16GB

Approximate run time: 0.5-5 minutes per sample (depending on number of reads, length of reference, and available compute).

ARM processor support: True




## Install and run

These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore Nextflow will need to be installed before attempting to run the workflow.

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Socker or Singularity is installed. This is controlled by the [`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles) parameter as exemplified below.

It is not required to clone or download the git repository in order to run the workflow.
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This will pull the repository in to the assets folder of Nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-installation-qualification –help
```

A demo dataset is provided for testing of the workflow. It can be downloaded using:

```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-installation-qualification/wf-installation-qualification-demo.tar.gz
tar -xzvf wf-installation-qualification-demo.tar.gz
```

The workflow can be run with the demo data using:

```
nextflow run epi2me-labs/wf-installation-qualification \
    --fastq wf-installation-qualification-demo/fastq \
    --references wf-installation-qualification-demo/references \
    --sample_sheet wf-installation-qualification-demo/sample_sheet.csv \
    -profile standard
```

For further information about running a workflow on the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).




## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts either FASTQ or BAM files as input.

The FASTQ or BAM input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ or BAM file; (ii) the path to a top-level directory containing FASTQ or BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ or BAM files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| references | string | Path to a directory containing FASTA reference files. | Accepted file extensions are '.fasta', '.fna', '.ffn', '.faa', '.frn', '.fa', '.txt', '.fa.gz', '.fna.gz', '.frn.gz', '.ffn.gz', '.fasta.gz'. |  |
| tests_config | string | Path to a JSON file containing parameters for the qualification tests. Each barcode must be recorded in this file along with its target reference sequence and whether the sample is a negative_control. | See data/tests_config.json for details. |  |
| counts | string | Path to a CSV file containing expected counts as a control. | The expected counts CSV file must contain columns named 'reference' and 'expected_counts' in order to be valid. the 'reference' column should contain names matching the names of reference sequences within the fasta files provided using --references. |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |
| prefix | string | Optional prefix attached to each of the output filenames. | Output filename format will be `<prefix>-filename.ext`. |  |


### Advanced options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| depth_coverage | boolean | Calculate depth coverage statistics and include them in the report. | This step can be a computational bottleneck. Set this to false if your reference sequences are >50mb to speed things up. | True |
| minimap_preset | string | Pre-defined parameter sets for `minimap2`, covering most common use cases. | Available parameter sets are: 'dna' (`-ax map-ont`), 'rna' (`-ax splice -uf`). | dna |
| minimap_args | string | String of command line arguments to be passed on to `minimap2`. | This overrides the options defined by `--minimap_preset` and allows for running the alignment step in a more customized way. |  |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Number of CPU threads to use for the alignment step. | The alignment process will run with this many threads (note that the memory used by minimap2 scales with the number of threads). The total CPU resources used by the workflow are constrained by the executor configuration and can be modified by changing `nextflow.config` or supplying an additional config file. | 4 |
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | wf-installation-qualification-report.html | Report for all samples | aggregated |
| Combined references | combined-refs.fasta | FASTA file containing all input references. | aggregated |
| Combined references index | combined-refs.fasta.fai | Index file for combined references FASTA. | aggregated |
| Per-read alignment stats | {{ alias }}.readstats.tsv | Bamstats per-read output TSV file. | per-sample |
| Per-reference alignment stats | {{ alias }}.flagstat.tsv | Bamstats flagstat output TSV file. | per-sample |
| Alignments BAM file | {{ alias }}.sorted.aligned.bam | BAM file with alignments of filtered input reads against the combined references. | per-sample |
| Alignments index file | {{ alias }}.sorted.aligned.bam.bai | Index for alignments BAM file. | per-sample |
| Workflow checkpoints | checkpoints.json | Structured workflow checkpoints for internal/onward use. | aggregated |
| Workflow results | results.json | Structured workflow results for internal/onward use. | aggregated |




## Pipeline overview

### 1. Combine reference files

All reference files in the directory passed to `--references` are concatenated.

### 2. Align reads

Input reads are aligned against the combined reference with [Minimap2](https://github.com/lh3/minimap2). If BAM files are used as input (with `--bam`), only reads in files without a reference in the SAM header are aligned. For other BAM files this step is skipped.

### 3. Create alignment statistics

[Bamstats](https://github.com/epi2me-labs/fastcat#bamstats) is used to create per-read and per-reference alignment statistics from the BAM files.

### 4. Calculate depth of coverage

Depth of coverage along the reference sequences is determined with [Mosdepth](https://github.com/brentp/mosdepth) (using 200 windows per reference sequence). To speed up the workflow, this step can be skipped by adding `--depth-coverage false`.

### 5. Run qualification tests to compare statisticss to specified criteria

Criteria and thresholds listed in a tests_config JSON file are use to determine which samples have passed or failed the qualification checks. 
Test samples are required to be above the stated alignment thresholds for a specified reference sequence, whereas control samples are required to be below set limits. The type of sample that each barcode relates to (`test` or `negative_control`) is taken from the tests_cofig JSON file.




## Troubleshooting

+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ Please see [here](https://labs.epi2me.io/trouble-shooting/) for how to resolve some common Nextflow issues and [here](https://labs.epi2me.io/how-to-exits/) for how to interpret command exit codes.




## FAQ's

*I cannot select a single reference file in the EPI2ME desktop app.* - When running the workflow via the desktop app, you need to provide a directory with reference files. If you only have a single file, you can create a directory to place your reference file inside and select this with the reference input option.

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-installation-qualification/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts

- [How to align your data](https://labs.epi2me.io/how-to-align/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.




