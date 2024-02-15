#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2
nextflow.preview.recursion=true

include { fastq_ingress; xam_ingress; } from "./lib/ingress"
include { process_references } from "./subworkflows/process_references"
include {
    getParams;
} from './lib/common'


OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")
MINIMAP_ARGS_PRESETS = [
    "dna": "-ax map-ont -y",
    "rna": "-ax splice -uf -y"
]


///////////////////////
// Workflow checkpoints
///////////////////////

process accumulateCheckpoints {
    label "wf_common"
    cpus 1
    input:
        path data
        val metadata
        path definitions
    output:
        path "checkpoints_${task.index}.json"
        val metadata
        path definitions
    publishDir params.out_dir, mode: 'copy', overwrite: true, pattern: "checkpoints_${task.index}.json", saveAs: { 'checkpoints.json' }
    script:
        // If the data list is lenth 1 then nextflow makes it not a list
        def data = data instanceof List ? data: [data]
        // Set-up the output file
        output = "checkpoints_${task.index}.json"
        // The run metadata (sample sheet)
        metaJson = new JsonBuilder(metadata).toPrettyString()
        // The checkpoint data
        checkpoint_data = data.getAt(0)
        // Our 1st checkpoint will not have a checkpoints file created
        // and the length of the data array wil be 1. Any subsequent checkpoint
        // will have the data followed by all of the previous checkpoint files
        // we need the last.
        if (data.size() > 1) {
            checkpoints_file = "--checkpoints_file ${data.getAt(-1)}"
        } else {
            checkpoints_file = ""
        }
    """
    echo '${metaJson}' > metadata.json
    accumulate_checkpoints.py ${output} \
        --output_definitions ${definitions} \
        --checkpoint_data ${checkpoint_data} \
        --metadata metadata.json \
        ${checkpoints_file}
    """
}

process sample_preparationCheckpoint {
  label "wf_common"
  cpus 1
  memory "2 GB"
  input: 
    tuple val(meta), path(sample), val(stats)
  output:
    path "checkpoint_info.json", emit: checkpoint
  script:
    String status = 'complete'

    // make our checkpoint data
    def checkpoint_data = [[
            sample: "${meta.alias}",
            checkpoint_name: "sample_preparation",
            status: "complete"
        ]]

    checkJson = new JsonBuilder(checkpoint_data).toPrettyString()
    """
    echo '$checkJson' > checkpoint_info.json
    """
}

process alignmentCheckpoint {
  label "wf_common"
  cpus 1
  memory "2 GB"
  input:
    tuple val(meta), path(bam), path(bai)
  output:
    path "checkpoint_info.json", emit: checkpoint
  script:
    String status = 'complete'

    // make our checkpoint data
    def checkpoint_data = [[
            sample: "${meta.alias}",
            checkpoint_name: "alignment",
            status: "complete",
            files: [ "alignment": "./${bam}", "alignment-index": "./${bai}" ]
        ]]

    checkJson = new JsonBuilder(checkpoint_data).toPrettyString()
    """
    echo '$checkJson' > checkpoint_info.json
    """
}

process reportingCheckpoint {
  label "wf_common"
  cpus 1
  memory "2 GB"
  input:
    path report
  output:
    path "checkpoint_info.json", emit: checkpoint
  script:
    String status = 'complete'

    // make our checkpoint data
    def checkpoint_data = [[
            sample: "",
            checkpoint_name: "reporting",
            status: "complete",
            files: [ 
              "workflow-report": "./${report}",
              "workflow-results": "./results.json"]
        ]]

    checkJson = new JsonBuilder(checkpoint_data).toPrettyString()
    """
    echo '$checkJson' > checkpoint_info.json
    """
}


///////////////////////
// Workflow processes
///////////////////////

process alignReads {
    label "wfalignment"
    cpus params.threads
    memory "12 GB"
    input:
        tuple val(meta), path(input)
        path combined_refs
        val is_xam
        val minimap_args
    output:
        tuple val(meta), path(bam_name)
    script:
        def sample_name = meta["alias"]
        bam_name = "${sample_name}.sorted.aligned.bam"
        int sorting_threads = Math.min((task.cpus / 3) as int, 3)
        int mapping_threads = task.cpus - sorting_threads
        // the minimum for `params.threads` in the schema is `4` and we should have
        // positive values for both thread vars, but can't hurt to make extra sure
        sorting_threads = Math.max(1, sorting_threads)
        mapping_threads = Math.max(1, mapping_threads)
    """
    ${is_xam ? "samtools fastq -T '*' $input" : "cat $input"} \
    | minimap2 -t $mapping_threads $minimap_args $combined_refs - \
    | samtools sort -@ ${sorting_threads - 1} -o $bam_name -
    """
}

process indexBam {
    label "wfalignment"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path(bam), path("*.bai")
    script:
    """
    samtools index $bam
    """
}

process bamstats {
    label "wfalignment"
    cpus 2
    memory "4 GB"
    input:
        tuple val(meta), path(bam), path(index)
    output:
        path "*.readstats.tsv", emit: read_stats
        path "*.flagstat.tsv", emit: flagstat
    script:
        def sample_name = meta["alias"]
    """
    bamstats $bam -s $sample_name -u -f ${sample_name}.flagstat.tsv -t $task.cpus \
    > ${sample_name}.readstats.tsv
    """
}

process addStepsColumn {
    // TODO: we don't need 200 windows for very short references; find heuristics for
    // determining window length / number for such cases
    label "wfalignment"
    cpus 1
    memory "2 GB"
    input: path "lengths.tsv"
    output: path "lengths_with_steps.tsv"
    """
    #!/usr/bin/env python
    import pandas as pd
    all = pd.read_csv('lengths.tsv', sep='\\t')
    all["step"] = all["lengths"]//200
    all = all.replace(0, 1)
    all.to_csv('lengths_with_steps.tsv', index=False, header=False, sep='\\t')
    """
}

process readDepthPerRef {
    // TODO: check if parallelisation with `xargs` or `parallel` is more efficient
    label "wfalignment"
    cpus 3
    memory "8 GB"
    input:
        tuple val(meta), path(alignment), path(index)
        path ref_len
    output:
        path outfname
    script:
        def sample_name = meta["alias"]
        outfname = "${sample_name}.all_regions.bed.gz"
    """
    while IFS=\$'\\t' read -r name lengths steps; do
        mosdepth -n --fast-mode --by "\$steps" --chrom "\$name" -t $task.cpus \
            ${sample_name}."\$name".temp $alignment \
        || echo "No alignments for "\$name""
        [[ -f ${sample_name}."\$name".temp.regions.bed.gz ]] && \
            cat ${sample_name}."\$name".temp.regions.bed.gz >> $outfname
    done < $ref_len

    # remove all the temp files
    find -name '${sample_name}.*.temp*' -delete
    """
}

process makeReport {
    label "wfalignment"
    cpus 1
    memory "12 GB"
    input:
        path "readstats/*"
        path "flagstat/*"
        path "refnames/*"
        path depths, stageAs: "depths/*"
        path counts
        path versions
        path params
    output:
        path "*.html"
    script:
    String depth_args = "--depths_dir depths"
    // we need to check against `.baseName` here because Nextflow includes the staging
    // directory in the `.name` of a `TaskPath`
    if (!(depths instanceof List) && depths.baseName == OPTIONAL_FILE.name) {
        depth_args = ""
    }
    String counts_args = (counts.name == OPTIONAL_FILE.name) ? "" : "--counts $counts"
    """
    workflow-glue report \
        --name wf-installation-qualification \
        --stats_dir readstats \
        --flagstat_dir flagstat \
        --refnames_dir refnames \
        --versions $versions \
        --params $params \
        $depth_args \
        $counts_args
    """
}


process getVersions {
    label "wfalignment"
    cpus 1
    memory "2 GB"
    output:
        path "versions.txt"
    script:
    """
    python --version | tr -s ' ' ',' | tr '[:upper:]' '[:lower:]' > versions.txt
    seqkit version | sed 's/ /,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | (head -n 1 && exit 0) | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    ezcharts --version | sed 's/ /,/' >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    bgzip --version | head -n1 | sed -E 's/(.*) /\\1,/' >> versions.txt
    """
}


// workflow module
workflow pipeline {
    take:
        sample_data
        refs
        counts
        depth_coverage
        tests_config
    main:
        // get params & versions
        workflow_params = getParams()
        software_versions = getVersions()

        // handle references
        refs = process_references(params.references)

        sample_data = sample_data
        | map { meta, path, stats -> [meta, path] }
        String minimap_args

        if (params.bam) {
            ch_branched = sample_data.branch { meta, bam ->
                to_align: meta["is_unaligned"]
                aligned: true
            }
            ch_to_align = ch_branched.to_align
            // `xam_ingress` sorts the BAMs, so we don't have to
            bam = ch_branched.aligned
        } else {
            // FASTQ input
            ch_to_align = sample_data
            bam = Channel.empty()
        }

        // run minimap
        if (! MINIMAP_ARGS_PRESETS.containsKey(params.minimap_preset)) {
            error "'--minimap_preset' needs to be one of " +
                "${MINIMAP_ARGS_PRESETS.keySet()}."
        }
        minimap_args = params.minimap_args ?: \
            MINIMAP_ARGS_PRESETS[params.minimap_preset]
        bam = bam
        | mix(
            alignReads(ch_to_align, refs.combined, params.bam as boolean, minimap_args)
        )
        | indexBam

        alignment_checkpoint = alignmentCheckpoint(bam)

        // get stats
        stats = bamstats(bam)

        // determine read_depth per reference / bam file if requested
        depth_per_ref = Channel.of(OPTIONAL_FILE)
        if (depth_coverage) {
            // add step column to ref lengths
            ref_lengths_with_steps = addStepsColumn(refs.lengths_combined)
            depth_per_ref = readDepthPerRef(bam, ref_lengths_with_steps)
        }
        read_stats = stats.read_stats.collect()
        flagstat = stats.flagstat.collect()
        
        qualification_tests = process_tests(
            read_stats,
            flagstat,
            refs.names_per_ref_file.collect(),
            refs.lengths_per_ref_file.collect(),
            tests_config
        )

        metadata = sample_data.map { meta, reads -> meta }.toList()

        // process the results into a json file for dissemination
        collected_results = collect_results(
            read_stats,
            qualification_tests,
            metadata)
        
        report = makeReport(
            read_stats,
            flagstat,
            refs.names_per_ref_file.collect(),
            depth_per_ref.collect(),
            counts,
            software_versions,
            workflow_params,
        )
    // If the emit values are altered, check that the output() process 
    // is run only on the appropriate outputs. 
    // Eg. output(jb2_conf.concat(*results[0..-3]))
    emit:
        alignments = bam.map { it[1] }
        indices = bam.map{ it[2] }
        per_read_stats = stats.read_stats
        per_file_stats = stats.flagstat
        report
        params_json = workflow_params
        software_versions
        combined_ref = refs.combined
        combined_ref_index = refs.combined_index
        collected_results = collected_results
        alignment_checkpoint = alignment_checkpoint
        metadata = metadata
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    label "wfalignment"
    cpus 1
    memory "12 GB"
    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*", saveAs: {
        f -> params.prefix ? "${params.prefix}-${f}" : "${f}" }
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}


process configure_jbrowse {
    label "wfalignment"
    cpus 1
    memory "2 GB"
    input:
        path(alignments)
        path(indexes)
        path(reference)
        path(ref_idx)
    output:
        path("jbrowse.json")
    script:
    ArrayList alignment_args = []
    int i = 0;
    for(a in alignments) {
        // don't be fooled into iterating over bam.size() here
        // when the cardinality is 1, bam.size() returns the filesize of the bam!
        this_bam = a
        this_bai = indexes[i]
        alignment_args << "--alignment '${params.out_dir}/${this_bam.name}' '${params.out_dir}/${this_bai.name}'"
        i++;
    }
    String alignment_args_str = alignment_args.join(' ')
    """
    workflow-glue configure_jbrowse \
        --reference '${reference}' '${params.out_dir}/${reference.name}' '${params.out_dir}/${ref_idx.name}' \
        ${alignment_args_str} > jbrowse.json
    """
}


process process_tests {
    label "wfalignment"
    cpus 1
    memory "2 GB"
    input:
        path "readstats/*"
        path "flagstat/*"
        path "refnames/*"
        path "reflengths/*"
        path "tests_config"
    output:
        path "qualification_tests.json"
    script:
    """
    workflow-glue process_tests \
        --stats_dir readstats \
        --flagstat_dir flagstat \
        --refnames_dir refnames \
        --reflengths_dir reflengths \
        --tests_config tests_config
    """
    
}


process collect_results {
    label "wfalignment"
    cpus 1
    memory "12 GB"
    input:
        path "readstats/*"
        path "qualification_tests"
        val "metadata"
    output:
        path "results.json"
    script:
    metaJson = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metaJson}' > metadata.json
    workflow-glue collect_results \
        --output "results.json" \
        --stats_dir readstats \
        --qualification_tests $qualification_tests \
        --meta "metadata.json"
    """
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)

    File checkpoints_file = new File("checkpoints.json");  

    if (checkpoints_file.exists() == true && workflow.resume == false){
        checkpoints_file.delete()
    } 

    // set tests_config file
    // (can't reassign a params value so create a local variable)
    if (!params.tests_config){
        tests_config = projectDir.resolve("./data/tests_config.json").toString()
    } else {
        tests_config = params.tests_config
    }



    Map ingress_args = [
        "sample_sheet": params.sample_sheet,
        "analyse_unclassified": params.analyse_unclassified,
        "stats": false,
    ]

    // get input data
    if (params.fastq) {
        sample_data = fastq_ingress(ingress_args + ["input": params.fastq])
    } else {
        sample_data = xam_ingress(
            ingress_args + ["input": params.bam, "keep_unaligned": true]
        )
    }

    counts = file(params.counts ?: OPTIONAL_FILE, checkIfExists: true)
    sampleprep = sample_preparationCheckpoint(sample_data)

    // Run pipeline
    results = pipeline(
        sample_data, params.references, counts, params.depth_coverage, tests_config
    )

    // create jbrowse file
    jb2_conf = configure_jbrowse(
        results.alignments.collect(),
        results.indices.collect(),
        results.combined_ref,
        results.combined_ref_index
    )


    // output all but the last two items of results
    output(jb2_conf.concat(*results[0..-3]))

    // The output definition file
    definitions = projectDir.resolve("./output_definition.json").toString()

    // checkpoints
    reporting_checkpoint = reportingCheckpoint(results.report)
    accumulateCheckpoints.scan(
      sampleprep.mix(*[results.alignment_checkpoint, reporting_checkpoint]
      ), results.metadata, definitions) 

}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
