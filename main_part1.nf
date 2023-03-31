#!/usr/bin/env nextflow

nextflow.enable.dsl=2

fastq_source_folder = Channel.fromPath(params.fastq_source_folder, type: "dir", checkIfExists: true).toList()
fq_list = file(params.fq_list)
reference_folder = file(params.reference_folder)
bed_file = file(params.bed_file)

knownsite1vfc = file(params.knownsite1vcf)
knownsite1vfctbi = file(params.knownsite1vcftbi)

knownsite2vfc = file(params.knownsite2vcf)
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

    pbrun bammetrics \
        --ref ${reference_folder}/${reference_name} \
        --bam ${output}/mark_dups_HG002.bam \
        --minimum-base-quality 2 \
        --minimum-mapping-quality 5 \
        --coverage-cap 250 \
        --interval-file ${bed_file} \
        --out-metrics-file ${output}/wgsmetrics_HConly.txt \
        --num-threads 32 \
        --logfile ${output}/bammetrics.log
        
    mkdir ${output}/collectmultiplemetrics -p

    pbrun collectmultiplemetrics \
        --ref ${reference_folder}/${reference_name} \
        --bam ${output}/mark_dups_HG002.bam \
        --out-qc-metrics-dir ${output}/collectmultiplemetrics \
        --gen-all-metrics \
        --bam-decompressor-threads 10 \
        --logfile ${output}/collectmultiplemetrics.log

    pbrun bqsr \
        --ref ${reference_folder}/${reference_name} \
        --in-bam ${output}/mark_dups_HG002.bam \
        --knownSites ${knownsite2vfc} \
        --knownSites ${knownsite1vfc} \
        --interval-file ${bed_file} \
        --out-recal-file ${output}/bqsr_recalibration_table.txt \
        --logfile ${output}/bqsr.log


    """
}
