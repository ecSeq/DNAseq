#!/usr/bin/env nextflow
// This file defines individual processes (separated for portability)

// perform trimming with cutadapt
process "cutadapt" {

    label "low"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    publishDir "${params.output}", pattern: "clipped/*.fastq.gz", mode: 'copy'
    publishDir "${params.output}", pattern: "clipped/logs/*.log", mode: 'move'

    input:
    tuple val(sample), path(reads)
    // eg. [sample, [read1.fastq.gz, read2.fastq.gz]]
    // eg. [sample, reads.fastq.gz]

    output:
    tuple val(sample), path("clipped/*.fastq.gz")
    // eg. [sample, [/path/to/clipped/sample_1.fastq, /path/to/clipped/sample_2.fastq]]
    // eg. [sample, /path/to/clipped/sample.fastq]
    path "clipped/*.fastq.gz"
    path "clipped/logs/*.log"

    script:
    if( params.SE )
        """
        mkdir clipped clipped/logs
        cutadapt -j ${task.cpus} -a ${params.forward} \\
        -q ${params.minQual} -m ${params.minLeng} -O ${params.minOver} \\
        -o clipped/${reads} ${reads} \\
        > clipped/logs/${sample}.log 2>&1
        """
    else
        """
        mkdir clipped clipped/logs
        cutadapt -j ${task.cpus} -a ${params.forward} -A ${params.reverse} \\
        -q ${params.minQual} -m ${params.minLeng} -O ${params.minOver} \\
        -o clipped/${reads[0]} \\
        -p clipped/${reads[1]} ${reads} \\
        > clipped/logs/${sample}.log 2>&1
        """
}


// generate quality control report
process "FastQC" {

    label "low"
    label "ignore"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    publishDir "${params.output}", pattern: "clipped/fastqc/*.{html,zip}", mode: 'move'
    publishDir "${params.output}", pattern: "clipped/fastqc/logs/*.log", mode: 'move'

    input:
    tuple val(sample), path(reads)
    // eg. [sample, [read1.fastq.gz, read2.fastq.gz]]
    // eg. [sample, reads.fastq.gz]

    output:
    path "clipped/fastqc/*.{html,zip}"
    path "clipped/fastqc/logs/${sample}.log"

    script:
    """
    mkdir clipped clipped/fastqc clipped/fastqc/logs
    fastqc ${reads} -threads ${task.cpus} \\
    -outdir=clipped/fastqc > clipped/fastqc/logs/${sample}.log
    """
}



// index the genome
process "indexing" {

    label "low"
    label "finish"
    tag "$type"

    maxForks "${params.fork}".toInteger()

    input:
    val type
    // eg. [BWA]
    path fasta
    path fai

    output:
    tuple val(type), path("genome/${type}")
    // eg. [BWA, /path/to/genome/BWA]

    script:
    if( type == "segemehl" )
        """
        mkdir genome genome/segemehl
        segemehl.x -d ${fasta} -x genome/segemehl/${fasta.baseName}.idx
        """
    else if( type == "STAR" )
        """
        mkdir genome genome/STAR
        STAR --runMode genomeGenerate \\
        --genomeDir genome/STAR/ \\
        --outFileNamePrefix genome/STAR/${fasta.baseName}. \\
        --genomeFastaFiles ${fasta} --runThreadN ${task.cpus}
        """
    else if( type == "bowtie2" )
        """
        mkdir genome genome/bowtie2
        bowtie2-build ${fasta} genome/bowtie2/${fasta.baseName}
        """
    else if( type == "BWA" || type == "BWA_MEM" )
        """
        mkdir genome genome/${type}
        cp ${fasta} genome/${type}
        bwa index genome/${type}/${fasta}
        """
}



// align the reads to the genome
process "mapping" {

    label "high"
    label "finish"
    tag "$type"

    maxForks "${params.fork}".toInteger()

    publishDir "${params.output}", pattern: "mapping/*.bam", mode: 'move'
    publishDir "${params.output}", pattern: "mapping/logs/*.log", mode: 'move'

    input:
    tuple val(sample), path(reads), val(type), path(index)
    // eg. [sample, [/path/to/clipped/sample_1.fastq, /path/to/clipped/sample_2.fastq], BWA, /path/to/genome/BWA]
    // eg. [sample, /path/to/clipped/sample.fastq, BWA, /path/to/genome/BWA]
    path fasta
    path fai

    output:
    path "mapping/${sample}.${type}.bam"
    path "mapping/logs/${sample}.${type}.log"

    script:
    if( type == "segemehl" )
        """
        mkdir mapping mapping/logs
        segemehl.x -d ${fasta} -i ${index}/${fasta.baseName}.idx \\
        ${params.SE ? "-q ${reads} " : "-q ${reads[0]} -p ${reads[1]} "} -t ${task.cpus} 2> mapping/logs/${sample}.${type}.log |
        samtools sort -O BAM -T segemehl.deleteme > mapping/${sample}.${type}.bam 
        """
    else if( type == "STAR" )
        """
        mkdir mapping mapping/logs
        STAR --genomeDir ${index}/ --readFilesIn ${reads} \\
        --outFileNamePrefix ${sample}.${type}. --outSAMattributes NH HI AS nM NM MD jM jI MC \\
        --runThreadN ${task.cpus} --readFilesCommand zcat --outStd SAM 2> mapping/logs/${sample}.${type}.log |
        samtools sort -O BAM -T STAR.deleteme > mapping/${sample}.${type}.bam
        """
    else if( type == "bowtie2" )
        """
        mkdir mapping mapping/logs
        bowtie2 -x ${index}/${fasta.baseName} \\
        ${params.SE ? "-U ${reads} " : "-1 ${reads[0]} -2 ${reads[1]} "} -p ${task.cpus} 2> mapping/logs/${sample}.${type}.log |
        samtools sort -O BAM -T bowtie2.deleteme > mapping/${sample}.${type}.bam 
        """
    else if( type == "BWA" )
        """
        mkdir mapping mapping/logs
        echo ${reads} | tr " " "\\n" |
        xargs -n2 -P2 -i sh -c "bwa aln ${index}/${fasta} '{}' > \$(basename '{}' .fastq).sai" || exit \$?
        bwa sampe ${index}/${fasta} *.sai ${reads} 2> mapping/logs/${sample}.${type}.log | 
        samtools sort -O BAM -T BWA.deleteme > mapping/${sample}.${type}.bam
        """
    else if( type == "BWA_MEM" )
        """
        mkdir mapping mapping/logs
        bwa mem ${index}/${fasta} ${reads} 2> mapping/logs/${sample}.${type}.log |
        samtools sort -O BAM -T BWA_MEM.deleteme > mapping/${sample}.${type}.bam
        """
}
