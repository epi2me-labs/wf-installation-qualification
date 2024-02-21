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
