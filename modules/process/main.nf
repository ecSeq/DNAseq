// This file defines individual processes (separated for portability)

// perform trimming with cutadapt
process "cutadapt" {

    label "DNAseq"
    label "low"
    label "finish"

    tag "${sample}"
    publishDir "${params.output}", pattern: "trimming/*.fastq.gz", mode: 'copy', enabled: params.keepReads ? true : false
    publishDir "${params.output}", pattern: "trimming/logs/${sample}.cutadapt.log", mode: 'move'

    maxForks "${params.fork}".toInteger()

    input:
    tuple val(sample), path(reads)
    // eg. [sample, [read1.fastq.gz, read2.fastq.gz]]
    // eg. [sample, reads.fastq.gz]

    output:
    tuple val(sample), path("trimming/*.fastq{,.gz}")
    // eg. [sample, [/path/to/trimming/sample_1.fastq, /path/to/trimming/sample_2.fastq]]
    // eg. [sample, /path/to/trimming/sample.fastq]
    path "trimming/*.fastq{,.gz}"
    path "trimming/logs/${sample}.cutadapt.log"

    script:
    if( params.SE )
        """
        mkdir trimming trimming/logs
        cutadapt -j ${task.cpus} -a ${params.forward} \\
        -q ${params.minQual} -m ${params.minLeng} -O ${params.minOver} \\
        -o trimming/${reads} ${reads} \\
        > trimming/logs/${sample}.cutadapt.log 2>&1
        """
    else
        """
        mkdir trimming trimming/logs
        cutadapt -j ${task.cpus} -a ${params.forward} -A ${params.reverse} \\
        -q ${params.minQual} -m ${params.minLeng} -O ${params.minOver} \\
        -o trimming/${reads[0]} \\
        -p trimming/${reads[1]} ${reads} \\
        > trimming/logs/${sample}.cutadapt.log 2>&1
        """
}


// generate quality control report
process "FastQC" {

    label "DNAseq"
    label "low"
    label "ignore"

    tag "${sample}"
    publishDir "${params.output}", pattern: "trimming/*.{html,zip}", mode: 'move'
    publishDir "${params.output}", pattern: "trimming/logs/${sample}.fastqc.log", mode: 'move'

    maxForks "${params.fork}".toInteger()

    input:
    tuple val(sample), path(reads)
    // eg. [sample, [read1.fastq.gz, read2.fastq.gz]]
    // eg. [sample, reads.fastq.gz]

    output:
    path "trimming/*.{html,zip}"
    path "trimming/logs/${sample}.fastqc.log"

    when:
    params.FastQC

    script:
    """
    mkdir trimming trimming/logs
    fastqc ${reads} -threads ${task.cpus} \\
    -outdir=trimming > trimming/logs/${sample}.fastqc.log 2>&1
    """
}



// index the genome
process "bowtie2_index" {

    label "DNAseq"
    label "low"
    label "finish"

    maxForks "${params.fork}".toInteger()

    input:
    path fasta
    path fai

    output:
    path "genome/bowtie2"

    script:
    """
    mkdir genome genome/bowtie2
    bowtie2-build ${fasta} genome/bowtie2/${fasta.baseName}
    """
}



// align the reads to the genome
process "bowtie2" {

    label "DNAseq"
    label "low"
    label "finish"

    tag "${sample}"
    publishDir "${params.output}", pattern: "mapping/*.bam", mode: params.bamQC ? 'copy' : 'move'
    publishDir "${params.output}", pattern: "mapping/logs/*.log", mode: 'move'

    maxForks "${params.fork}".toInteger()

    input:
    tuple val(sample), path(reads), path(index)
    // eg. [sample, [/path/to/clipped/sample_1.fastq, /path/to/clipped/sample_2.fastq], /path/to/genome/BWA]
    // eg. [sample, /path/to/clipped/sample.fastq, /path/to/genome/BWA]
    path fasta
    path fai

    output:
    tuple val(sample), path("mapping/${sample}.bowtie2.bam")
    // eg. [sample, /path/to/mapping/sample.bowtie2.bam]
    path "mapping/logs/${sample}.bowtie2.log"

    script:
    """
    mkdir mapping mapping/logs
    bowtie2 -x ${index}/${fasta.baseName} \\
    ${params.SE ? "-U ${reads} " : "-1 ${reads[0]} -2 ${reads[1]} "} \\
    -p ${task.cpus} 2> mapping/logs/${sample}.bowtie2.log |
    samtools sort -O BAM -T bowtie2.deleteme > mapping/${sample}.bowtie2.bam 
    """
}


// perform QC of alignments
process "bamQC" {

    label "DNAseq"
    label "low"
    label "finish"

    tag "${sample}"
    publishDir "${params.output}", pattern: "mapping/${sample}/*.{txt,pdf}", mode: 'move'
    publishDir "${params.output}", pattern: "mapping/logs/${sample}.bamqc.log", mode: 'move'

    maxForks "${params.fork}".toInteger()

    input:
    tuple val(sample), path(bam)
    // eg. [sample, /path/to/mapping/sample.bowtie2.bam]

    output:
    path "mapping/${sample}/*.{txt,pdf}"
    path "mapping/logs/${sample}.bamqc.log"

    when:
    params.bamQC

    script:
    """
    mkdir mapping mapping/logs
    qualimap bamqc -bam ${bam} -outdir mapping/${sample} -outformat pdf > mapping/logs/${sample}.bamqc.log 2>&1
    """
}
