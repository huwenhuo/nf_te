#!/usr/bin/env nextflow

nextflow.enable.dsl=2

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_name, file(row.bam_ori_file)) }
    .set { sample_channel }

process sortBamByQueryName {

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("sorted_${sample_id}.bam")  // Output as tuple

    script:
    """
    singularity exec ${params.img_samtools} samtools sort -@ 9 -m 5G -n -o sorted_${sample_id}.bam ${bam_file} 
    sample_name=\$(samtools view -H ${bam_file} | grep SM | awk '{print \$3}' | sed 's/SM://')
    echo ${sample_id} \${sample_name} > id.txt
    mkdir -p ${params.output_dir}/${sample_id}
    cp id.txt ${params.output_dir}/${sample_id}/id.txt

    """
}


// Process for converting BAM to FASTQ using samtools fastq
process bamToFastq {

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq"), path("${sample_id}_R2.fastq")  // Output as tuple

    script:
    """
    singularity exec ${params.img_samtools} samtools fastq -@ 9 ${sorted_bam} -1 ${sample_id}_R1.fastq -2 ${sample_id}_R2.fastq
    """
}


// Define process for STAR alignment
process starAlignment {

    input:
    tuple val(sample_id), path(R1_fastq), path(R2_fastq)


    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    """
    singularity exec ${params.img_star} STAR --genomeDir ${params.genomeDir} \\
         --runThreadN 28 \\
         --outSAMtype BAM SortedByCoordinate \\
         --outFilterMultimapNmax 1000 \\
         --outSAMmultNmax 1 \\
         --outFilterMismatchNmax 3 \\
         --outMultimapperOrder Random \\
         --winAnchorMultimapNmax 1000 \\
         --alignEndsType EndToEnd \\
         --alignIntronMax 1 \\
         --alignMatesGapMax 350 \\
         --quantMode TranscriptomeSAM \\
         --readFilesIn ${R1_fastq} ${R2_fastq} \\
         --readFilesCommand cat \\
         --outFileNamePrefix ${sample_id}_
    """

}

process teCountProcess {

    // Input files (e.g., BAM files and GTF files)
    input:
    tuple val(sample_id), path(bam_file)
    path GTF_FILE
    path TE_GTF_FILE

    // Output (with sample_id included in the output directory)
    output:
    tuple val(sample_id), path("tecount.cntTable")

    // Publish directory (with sample_id)
    publishDir "${params.output_dir}/${sample_id}/", mode: 'copy'


    script:
    """
    singularity exec ${params.img_tecount} TEcount \\
            --sortByPos --format BAM --mode multi \\
            -b ${bam_file} \\
            --GTF ${GTF_FILE} \\
            --TE ${TE_GTF_FILE} \\
            --project tecount
    touch tecount.cntTable
    """
}

process teLocalProcess {

    input:
    tuple val(sample_id), path(bam_file)
    path GTF_FILE
    path TE_loc_GTF_FILE

    output:
    tuple val(sample_id), path("telocal.cntTable") 

    // Publish directory (with sample_id)
    publishDir "${params.output_dir}/${sample_id}/", mode: 'copy'


    script:
    """
    singularity exec ${params.img_telocal} TElocal --sortByPos -b ${bam_file} --GTF ${GTF_FILE} --TE ${TE_loc_GTF_FILE} --stranded reverse --project telocal
    touch telocal.cntTable
    """
}


// Workflow to run all the processes
workflow {

	sorted_bam = sortBamByQueryName(sample_channel)
	fastq_files = bamToFastq(sorted_bam)
	aligned_bam = starAlignment(fastq_files)
	teCountProcess(aligned_bam, params.GTF_FILE, params.TE_GTF_FILE)
	teLocalProcess(aligned_bam, params.GTF_FILE, params.TE_loc_GTF_FILE )

}

