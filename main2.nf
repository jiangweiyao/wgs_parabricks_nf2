#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//params.input_fastqs = "/home/ubuntu/parabricks_sample/Data/*_{1,2}.fq.gz"
//fastq_files = Channel.fromFilePairs(params.input_fastqs, type: 'file')

params.fastq_source_folder = "/home/ubuntu/WGS_V2/G001_0056"
fastq_source_folder = Channel.fromPath(params.fastq_source_folder, type: "dir", checkIfExists: true).toList()

params.fq_list = "/home/ubuntu/WGS_V2/fastq_list.txt"
fq_list = file(params.fq_list)

params.reference_folder = "/home/ubuntu/WGS_V2/ucsc_goldenPath_analysisSet_grch38_HG002_altrevised_reference_v4.2.1"
reference_folder = file(params.reference_folder)
params.reference_name = "HG002_hg38.analysisSet.fa"
params.out = "parabricks_output2_patched"
params.output = "output_patched"

params.bed_file = "/home/ubuntu/WGS_V2/GIAB_benchmark_files/v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
bed_file = file(params.bed_file)

params.knownsite1vcf = "/home/ubuntu/WGS_V2/GIAB_benchmark_files/v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
knownsite1vfc = file(params.knownsite1vcf)
params.knownsite1vcftbi = "/home/ubuntu/WGS_V2/GIAB_benchmark_files/v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
knownsite1vfctbi = file(params.knownsite1vcftbi)

params.knownsite2vcf = "/home/ubuntu/WGS_V2/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
knownsite2vfc = file(params.knownsite2vcf)
params.knownsite2vcftbi = "/home/ubuntu/WGS_V2/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
knownsite2vfctbi = file(params.knownsite2vcftbi)

workflow {
    fq2bam(fq_list, fastq_source_folder, reference_folder, params.reference_name, params.output, bed_file, knownsite1vfc, knownsite1vfctbi, knownsite2vfc, knownsite2vfctbi)
}

process fq2bam {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true
    cpus 32
    memory '120 GB'
    input:
    file(fq_list)
    file(fastq_source_folder)
    file(reference_folder)
    val(reference_name)
    val(output)
    file(bed_file)
    file(knownsite1vfc)
    file(knownsite1vfctbi)
    file(knownsite2vfc)
    file(knownsite2vfctbi)

    output:
    file("${output}/*")

    """
    mkdir ${output}

    pbrun fq2bam \
        --ref ${reference_folder}/${reference_name} \
        --in-fq-list ${fq_list} \
        --bwa-options="-M -Y" \
        --fix-mate \
        --out-duplicate-metrics ${output}/dupemetrics.txt \
        --out-qc-metrics-dir ${output}/qc_metrics \
        --optical-duplicate-pixel-distance 500 \
        --out-bam ${output}/mark_dups_HG002.bam \
        --logfile ${output}/fq2bam.log

    """
}
