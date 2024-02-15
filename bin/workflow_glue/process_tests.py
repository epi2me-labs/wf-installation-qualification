#!/usr/bin/env python

"""Process installation qualification test."""

import json

from .util import get_named_logger, wf_parser  # noqa: ABS101
from .report_utils import read_data  # noqa: ABS101


def main(args):
    """Run entry point."""
    logger = get_named_logger("Process tests")
    # read in files
    stats_df = read_data.bamstats(args.stats_dir)
    flagstat_df = read_data.flagstat(args.flagstat_dir)
    tests_config = json.load(open(args.tests_config, 'rb'))
    refname2reffile = read_data.refnames(args.refnames_dir)
    reflengths = read_data.reflengths(args.reflengths_dir)

    # add a column with the respective ref. files to the stats dataframes
    for df in (stats_df, flagstat_df):
        if df is None:
            continue
        try:
            df["ref_file"] = (
                df["ref"]
                .apply(lambda ref: refname2reffile[ref])
                .astype(read_data.CATEGORICAL)
            )
        except KeyError as e:
            (missing_ref,) = e.args
            raise ValueError(
                f"Reference '{missing_ref}' not found in the provided "
                f"reference files."
            )

    # get tests results
    tests_out = {}
    for sample_name in list(tests_config['correct_references'].keys()) + \
        tests_config['control_samples']:  # noqa: E125
        found_fastq = sample_name in stats_df['sample_name'].values
        n_reads = stats_df.query(f"sample_name == '{sample_name}'").shape[0]
        n_bases = stats_df.query(f"sample_name == '{sample_name}'")['read_length'].sum()
        # get total alignments
        sample_flagstat_df = flagstat_df.query(f"sample_name == '{sample_name}'")
        alignments = (
            sample_flagstat_df.groupby("ref_file")["primary"].sum().drop("unmapped")
            ).sum()
        perc_alignments = alignments / n_reads * 100
        tests_out[sample_name] = {}
        tests_out[sample_name]['found_fastq'] = found_fastq
        tests_out[sample_name]['n_reads'] = n_reads
        tests_out[sample_name]['n_bases'] = int(n_bases)
        tests_out[sample_name]['perc_alignments'] = perc_alignments

    for sample_name in tests_config['correct_references'].keys():
        # get alignments per ref file (the table will have one row per ref file)
        sample_flagstat_df = flagstat_df.query(f"sample_name == '{sample_name}'")
        alignments_per_ref_file = (
            sample_flagstat_df.groupby("ref_file")["primary"].sum().drop("unmapped")
            ).sort_index()
        correct_ref = tests_config['correct_references'][sample_name]
        ref_found = False
        for ref in alignments_per_ref_file.index:
            if correct_ref in ref:
                perc_correct = alignments_per_ref_file[ref] / \
                    tests_out[sample_name]['n_reads'] * 100
                ref_length = reflengths[correct_ref]
                coverage = tests_out[sample_name]['n_bases'] / ref_length
                ref_found = True
        if ref_found is False:
            raise ValueError(
                f"Correct reference '{correct_ref}*' not found in the provided "
                f"reference files."
                )
        tests_out[sample_name]['perc_correct'] = perc_correct
        tests_out[sample_name]['coverage'] = round(coverage, 2)

        # test against thresholds
        if tests_out[sample_name]['n_reads'] >= \
            tests_config['parameters']["read_count_threshold"]:  # noqa: E125
            tests_out[sample_name]['read_count_threshold'] = 'PASS'
        else:
            tests_out[sample_name]['read_count_threshold'] = 'FAIL'
        if perc_correct/100 >= tests_config['parameters'][
            "positive_correctly_mapped_threshold"]:  # noqa: E125
            tests_out[sample_name]['positive_correctly_mapped_threshold'] = 'PASS'
        else:
            tests_out[sample_name]['positive_correctly_mapped_threshold'] = 'FAIL'
        if coverage >= tests_config['parameters']["positive_coverage_threshold"]:
            tests_out[sample_name]['positive_coverage_threshold'] = 'PASS'
        else:
            tests_out[sample_name]['positive_coverage_threshold'] = 'FAIL'
        tests_out[sample_name]['control_coverage_limit'] = 'NA'
        tests_out[sample_name]['control_mapped_limit'] = 'NA'
        tests_out[sample_name]['control_coverage_limit'] = 'NA'

    control_ref_length = reflengths[tests_config['parameters']
        ["control_coverage_reference"]]  # noqa: E128
    for sample_name in tests_config['control_samples']:
        coverage = tests_out[sample_name]['n_bases'] / control_ref_length
        # test against thresholds
        if tests_out[sample_name]['n_reads'] > \
            tests_config['parameters']["control_mapped_limit"]:  # noqa: E125
            tests_out[sample_name]['control_mapped_limit'] = 'FAIL'
        else:
            tests_out[sample_name]['control_mapped_limit'] = 'PASS'
        if tests_out[sample_name]['perc_alignments']/100 > \
            tests_config['parameters']["control_mapped_limit"]:  # noqa: E125
            tests_out[sample_name]['control_mapped_limit'] = 'FAIL'
        else:
            tests_out[sample_name]['control_mapped_limit'] = 'PASS'
        if coverage <= tests_config['parameters']["control_coverage_limit"]:
            tests_out[sample_name]['control_coverage_limit'] = 'PASS'
        else:
            tests_out[sample_name]['control_coverage_limit'] = 'FAIL'
        tests_out[sample_name]['positive_coverage_threshold'] = 'NA'
        tests_out[sample_name]['perc_correct'] = 'NA'
        tests_out[sample_name]['read_count_threshold'] = 'NA'
        tests_out[sample_name]['positive_correctly_mapped_threshold'] = 'NA'
        tests_out[sample_name]['coverage'] = round(coverage, 2)

    for sample_name in tests_out:
        if 'FAIL' in [tests_out[sample_name][x] for x in tests_out[sample_name]]:
            tests_out[sample_name]['sample_pass'] = 'FAIL'
        else:
            tests_out[sample_name]['sample_pass'] = 'PASS'

    json.dump(tests_out, open('qualification_tests.json', 'w'), indent=4)
    logger.info("Qualification tests written to qualification_tests.json.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("tests")
    parser.add_argument(
        "--stats_dir",
        help="directory with `bamstats` per-read stats",
    )
    parser.add_argument(
        "--flagstat_dir",
        help="directory with `bamstats` per-file stats",
    )
    parser.add_argument(
        "--refnames_dir",
        help="directory with files containing reference names",
    )
    parser.add_argument(
        "--reflengths_dir",
        help="directory with files containing reference lengths",
    )
    parser.add_argument(
        "--tests_config",
        default=None,
        help="JSON file with quallification tests parameters",
    )
    return parser
