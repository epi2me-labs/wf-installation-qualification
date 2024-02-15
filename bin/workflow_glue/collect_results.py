#!/usr/bin/env python

"""Collect results into wf.Sample."""

import json

import workflow_glue.results_schema as wf
from .report_utils import read_data  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


def fastcat_stats(read_stats_df, alias):
    """Put fastcat stats into the model."""
    if alias in read_stats_df.groups.keys():
        alias_data = read_stats_df.get_group(alias)
        result = wf.FastqStats(
            n_seqs=len(alias_data['read_length']),
            n_bases=alias_data['read_length'].sum(),
            min_length=alias_data['read_length'].min(),
            max_length=alias_data['read_length'].max(),
            mean_quality=alias_data['mean_quality'].mean()
            )
    else:

        result = wf.FastqStats(
            n_seqs=0,
            n_bases=0,
            min_length=0,
            max_length=0,
            mean_quality=0
        )

    return result


def qual_tests(qualification_tests, alias):
    """Put qualification tests in model."""
    if alias in qualification_tests:
        sample_data = qualification_tests[alias]
        tests = wf.QualificationTests(
            found_fastq=sample_data['found_fastq'],
            n_reads=sample_data['n_reads'],
            perc_alignments=sample_data['perc_alignments'],
            coverage=sample_data['coverage'],
            read_count_threshold=sample_data['read_count_threshold'],
            positive_correctly_mapped_threshold=sample_data[
                'positive_correctly_mapped_threshold'],
            control_mapped_limit=sample_data['control_mapped_limit'],
            positive_coverage_threshold=sample_data['positive_coverage_threshold'],
            control_coverage_limit=sample_data['control_coverage_limit'],
            sample_pass=sample_data['sample_pass']
        )
    else:
        tests = wf.QualificationTests(
            sample_pass="FAIL",
        )
    return tests


def main(args):
    """Run entry point."""
    logger = get_named_logger("Collect results")

    qualification_tests = json.load(open(args.qualification_tests, 'rb'))
    read_stats_df = read_data.bamstats(args.stats_dir)
    read_stats_df['sample_name'] = read_stats_df['sample_name'].astype(str)
    read_stats_df = read_stats_df.groupby('sample_name')
    meta = json.load(open(args.meta, 'rb'))

    # collect results for each sample
    samples = []
    w_pass = True
    for sample in meta:
        s_pass = True
        alias = sample['alias']
        sample_type = sample['type']
        barcode = sample['barcode']
        fastcat = fastcat_stats(read_stats_df, alias)
        quals = qual_tests(qualification_tests, alias)

        # set w_pass to false if any samples fail
        if quals.sample_pass == 'FAIL':
            w_pass = False
            s_pass = False
        # create results and samples
        results = wf.ResultsContents(
            tests=quals,
            fastq=fastcat
        )

        # Record if sample passes
        sample_checks = wf.CheckResult(
            check_name='all_qualification_tests_pass',
            check_pass=s_pass
            )

        sample_mod = wf.Sample(
            alias=alias,
            sample_type=sample_type,
            sample_pass=(quals.sample_pass == 'PASS'),
            results=results,
            barcode=barcode,
            sample_checks=[sample_checks]
        )
        samples.append(sample_mod)

    # Record if workflow passes
    workflow_checks = wf.CheckResult(
        check_name='all_sample_tests_pass',
        check_pass=w_pass
    )

    workflow = wf.WorkflowResult(
        samples=samples,
        workflow_checks=[workflow_checks],
        workflow_pass=w_pass
    )

    with open(args.output, 'w') as f:
        f.write(json.dumps(workflow.dict(), indent=4))

    logger.info(f"Results collected and written to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("collect_results")
    parser.add_argument(
        "--output",
        required=True,
        help="Output file to write JSON to",
    )
    parser.add_argument(
        "--stats_dir",
        required=True,
        help="Directory containing bamstats",
    )
    parser.add_argument(
        "--qualification_tests",
        required=True,
        help="JSON file containing qualification tests results",
    )
    parser.add_argument(
        "--meta",
        required=True,
        help="Samples metadata",
    )
    return parser
