#!/bin/bash
set -exo pipefail

get-test_data-from-aws () {
    # get aws-cli
    curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
    unzip -q awscliv2.zip

    # get test data
    aws/dist/aws s3 cp --recursive --quiet \
        "$S3_TEST_DATA" \
        test_data_from_S3
}

input_path=$1
input_type=$2
wf_output_dir=$3
sample_sheet=$4

# `input_path`, `ipnut_type`, and `wf_output_dir` are required
if ! [[ $# -eq 3 || $# -eq 4 ]]; then
    echo "Provide 2 or 3 arguments!" >&2
    exit 1
fi

# `input_type` needs to be either "fastq" or "bam"
if [[ $input_type != "fastq" && $input_type != "bam" ]]; then
    echo "The second argument must be 'fastq' or 'bam'!"
    exit 1
fi

# get test data from s3 if required
if [[ $input_path =~ ^s3:// ]]; then
    get-test_data-from-aws
    input_path="$PWD/test_data_from_S3/${input_path#*test_data/}"
    [[ -n $sample_sheet ]] &&
        sample_sheet="$PWD/test_data_from_S3/${sample_sheet#*test_data/}"
fi

# add CWD if paths are relative
[[ ( $input_path != /* ) ]] && input_path="$PWD/$input_path"
[[ ( $wf_output_dir != /* ) ]] && wf_output_dir="$PWD/$wf_output_dir"
[[ ( -n $sample_sheet ) && ( $sample_sheet != /* ) ]] &&
    sample_sheet="$PWD/$sample_sheet"

# add flags to parameters (`$input_path` could contain a space; so need to use an array
# here)
input_path=("--input" "$input_path")
input_type="--type $input_type"
wf_output_dir="--wf-output-dir $wf_output_dir"
[[ -n $sample_sheet ]] && sample_sheet="--sample_sheet $sample_sheet"

# get container hash from config
img_hash=$(grep 'common_sha.\?=' nextflow.config | grep -oE '(mr[0-9]+_)?sha[0-9,a-f,A-F]+')

# run test
docker run -v "$PWD":"$PWD" \
    ontresearch/wf-common:"$img_hash" \
    python "$PWD/test/test_ingress.py" "${input_path[@]}" $input_type $wf_output_dir $sample_sheet
