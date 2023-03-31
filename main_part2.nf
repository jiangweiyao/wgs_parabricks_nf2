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

ref_bed = file(params.ref_bed)

workflow {
    fq2bam(fq_list, fastq_source_folder, reference_folder, params.reference_name, params.output, bed_file, knownsite1vfc, knownsite1vfctbi, knownsite2vfc, knownsite2vfctbi)
    metrics(fq2bam.out[1], reference_folder, params.reference_name, ref_bed)
}

process fq2bam {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true
    cpus 32
    memory '121 GB'
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
    file("${output}/*.bam")    

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
        --out-bam ${output}/mark_dups.bam \
        --logfile ${output}/fq2bam.log

    """
}


process metrics {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true
    cpus 32
    memory '121 GB'
    input:
    file(input_bam)
    file(reference_folder)
    val(reference_name)
    file(ref_bed)
    

    output:
    file("mismatch_outs/*")

    """
mkdir -p mismatch_outs/ehists
mkdir -p mismatch_outs/mhists
mkdir -p mismatch_outs/ihists

awk -v OFS='\t' {'print \$1,\$2'} ${reference_folder}/${reference_name}.fai > ./grch38_bedtools.txt
cp ${ref_bed} tmp.bed
bedtools intersect -a ${input_bam} -b tmp.bed -sorted > sample_HConly.bam

### loop through GC bins
for num in \$(seq 0.20 0.05 0.70); do

    GCmax=\${num}
    GCmin=\$(bc <<< "\${GCmax}-0.04")
    GCmaxDecimal=\$(echo \${GCmax} | cut -d '.' -f 2)
    if [[ \${GCmaxDecimal} == "00" ]]; then GCmaxDecimal="100"; fi

    # need to create sub-BAM (or else stats will be on complete BAM)
    outbam=sample.0\${GCmin}to\${GCmax}.bam

    reformat.sh \
        -Xmx120g \
        in=sample_HConly.bam \
        out=\${outbam} \
        mingc=\${GCmin} \
        maxgc=\${GCmax} \
        mappedonly=t \
        primaryonly=t \
        ref=${reference_folder}/${reference_name} \
        crashjunk=f \
        tossjunk=t

    # now get statistics
    reformat.sh \
        -Xmx120g \
        in=\${outbam} \
        mingc=\${GCmin} \
        maxgc=\${GCmax} \
        mappedonly=t \
        primaryonly=t \
        ref=${reference_folder}/${reference_name} \
        ehist=ehist.sample.\${GCmin}to\${GCmaxDecimal}.txt \
        mhist=mhist.sample.\${GCmin}to\${GCmaxDecimal}.txt \
        indelhist=ihist.sample.\${GCmin}to\${GCmaxDecimal}.txt

    cp ehist.*.txt mismatch_outs/ehists/
    cp mhist.*.txt mismatch_outs/mhists/
    cp ihist.*.txt mismatch_outs/ihists/

    # clear BAM to maintain scratch space
    rm \${outbam}

done
    """
}
