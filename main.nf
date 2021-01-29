#!/usr/bin/env nextflow

// DSL2 BRANCH
nextflow.preview.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\

         ===========================================
          E C S E Q - D N A s e q   P I P E L I N E
         ===========================================
         ~ version ${workflow.manifest.version}

         Usage: 
              nextflow run epidiverse/template [OPTIONS]...

         Options: GENERAL
              --input [path/to/input/dir]     [REQUIRED] Provide the directory containing fastq file(s) in "*{1,2}.fastq.gz" format

              --reference [path/to/ref.fa]    [REQUIRED] Provide the path to the reference genome in fasta format

              --output [STR]                  A string that can be given to name the output directory. [default: "."]

              --SE                            Indicate to the pipeline whether fastq files are SE reads in "*.fastq.gz" format. [default: off]


         Options: ALIGNMENT [REQUIRED]
              --bowtie2                       Perform alignments with "bowtie2". [default: off]

              --BWA                           Perform alignments with "BWA". [default: off]

              --BWA_MEM                       Perform alignments with "BWA_MEM". [default: off]

              --segemehl                      Perform alignments with "segemehl". [default: off]

              --STAR                          Perform alignments with "STAR". [default: off]

         Options: ADDITIONAL
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    


         Example: 
              nextflow run ecseq/dna \
              --input /path/to/input/dir \
              --reference /path/to/genome.fa \
              --bowtie2 --BWA --BWA_MEM --segemehl --STAR

    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// PRINT VERSION AND EXIT
if(params.version){
    println """\
         ===========================================
          E C S E Q - D N A s e q   P I P E L I N E
         ===========================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}


// DEFINE PATHS # these are strings which are used to define input Channels,
// but they are specified here as they may be referenced in LOGGING
fasta = file("${params.reference}", checkIfExists: true, glob: false)
fai = file("${params.reference}.fai", checkIfExists: true, glob: false)
reads_path = params.SE ? "${params.input}/*.fastq.gz" : "${params.input}/*{1,2}.fastq.gz"


// populate list of alignment software based on user input params
if (!params.bowtie2 && !params.BWA && !params.BWA_MEM && !params.segemehl && !params.STAR) {
    error "ERROR: please specify one or more software to align with!"
} else {
    mapperList = []
    if (params.bowtie2){mapperList.add("bowtie2")}
    if (params.BWA){mapperList.add("BWA")}
    if (params.BWA_MEM){mapperList.add("BWA_MEM")}
    if (params.segemehl){mapperList.add("segemehl")}
    if (params.STAR){mapperList.add("STAR")}
}


// PRINT STANDARD LOGGING INFO
log.info ""
log.info "         ==========================================="
log.info "          E C S E Q - D N A s e q   P I P E L I N E"
if(params.debug){
log.info "         (debug mode enabled)"
log.info "         ===========================================" }
else {
log.info "         ===========================================" }
log.info "         ~ version ${workflow.manifest.version}"
log.info ""
log.info "         input dir    : ${workflow.profile.tokenize(",").contains("test") ? "-" : "${reads_path}"}"
log.info "         reference    : ${params.reference}"
log.info "         output dir   : ${params.output}"
log.info "         mapper(s)    : ${mapperList}"
log.info "         mode         : ${params.SE ? "single-end" : "paired-end"}"
log.info ""
log.info "         ================================================"
log.info "         RUN NAME: ${workflow.runName}"
log.info ""



////////////////////
// STAGE CHANNELS //
////////////////////

/*
 *   Channels are where you define the input for the different
 *    processes which make up a pipeline. Channels indicate
 *    the flow of data, i.e. the "route" that a file will take.
 */

// STAGE BAM FILES FROM TEST PROFILE # this establishes the test data to use with -profile test
if ( workflow.profile.tokenize(",").contains("test") ){

        include { check_test_data } from './libs/functions.nf' params(readPaths: params.readPaths, singleEnd: params.SE)
        READS = check_test_data(params.readPaths, params.SE)

} else {

    // STAGE READS CHANNELS # this defines the normal input when test profile is not in use
    READS = Channel
        .fromFilePairs(reads_path, size: params.SE ? 1 : 2)
        .ifEmpty{ exit 1, "ERROR: cannot find valid read files in dir: ${params.input}\n \
        The pipeline will expect PE reads in compressed *{1,2}.${params.extension} format\n \
        unless you have specified the --SE parameter or a different extension using --extension"}
        .map{ tuple(it[0], it[1]) }
        .take(params.take.toInteger())

}

// STAGE ALIGNMENT TYPES CHANNEL # this defines which software to run based on user parameters
TYPES = Channel.fromList(mapperList)



////////////////////
// BEGIN PIPELINE //
////////////////////

/*
 *   Workflows are where you define how different processes link together. They
 *    may be modularised into "sub-workflows" which must be named eg. 'DNAseq'
 *    and there must always be one MAIN workflow to link them together, which
 *    is always unnamed.
 */

// INCLUDES # here you must give the relevant process files from the libs directory 
include {cutadapt;FastQC;indexing;mapping} from './libs/process.nf' params(params)

// SUB-WORKFLOWS
workflow 'DNAseq' {

    // take the initial Channels and paths
    take:
        READS
        TYPES
        fasta
        fai

    // here we define the structure of our workflow i.e. how the different processes lead into each other
    // eg. process(input1, input2, etc.)
    // eg. process.out[0], process.out[1], etc.
    // index numbers [0],[1],etc. refer to different outputs defined for processes in process.nf
    // ALWAYS PAY ATTENTION TO CARDINALITY!!
    main:
        // first we will perform trimming with cutadapt
        cutadapt(READS)
        // the trimmed reads are then analysed with FastQC
        FastQC(cutadapt.out[0])

        // we can also already start indexing as these processes do not depend on each other
        indexing(TYPES,fasta,fai)

        // now we have a series of trimmed reads from cutadapt eg. [sample_name, [clipped_1.fastq.gz,clipped_2.fastq.gz]]
        // we also have a series of indexed genomes for mapping eg. [BWA, /path/to/index/dir]
        // so we should combine these channels to make all possible combinations:
        reads_and_software = cutadapt.out[0].combine(indexing.out)
        // eg. [sample_name, [clipped_1.fastq.gz,clipped_2.fastq.gz], BWA, /path/to/index/dir]
        // eg. [sample_name, [clipped_1.fastq.gz,clipped_2.fastq.gz], bowtie2, /path/to/index/dir]

        // and now we can run the alignments!
        mapping(reads_and_software,fasta,fai)

}

// MAIN WORKFLOW 
workflow {

    // call sub-workflows eg. WORKFLOW(Channel1, Channel2, Channel3, etc.)
    main:
        DNAseq(READS, TYPES, fasta, fai)

}


//////////////////
// END PIPELINE //
//////////////////

// WORKFLOW TRACING # what to display when the pipeline finishes
// eg. with errors
workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

// eg. in general
workflow.onComplete {

    log.info ""
    log.info "         Pipeline execution summary"
    log.info "         ---------------------------"
    log.info "         Name         : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}"
    log.info "         Profile      : ${workflow.profile}"
    log.info "         Launch dir   : ${workflow.launchDir}"    
    log.info "         Work dir     : ${workflow.workDir} ${!params.debug && workflow.success ? "(cleared)" : "" }"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""

    // run a small clean-up script to remove "work" directory after successful completion 
    if (!params.debug && workflow.success) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}