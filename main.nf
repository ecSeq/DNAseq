// PRINT HELP AND EXIT
if(params.help){
    include { printHelp } from './lib/functions.nf'
    printHelp()
}

// PRINT VERSION AND EXIT
if(params.version){
    include { printVersion } from './lib/functions.nf'
    printVersion()
}


// DEFINE PATHS # these are strings which are used to define input Channels,
// but they are specified here as they may be referenced in LOGGING
fasta = file("${params.reference}", checkIfExists: true, glob: false)
fai = file("${params.reference}.fai", checkIfExists: true, glob: false)
reads_path = params.SE ? "${params.input}/*.fastq" : "${params.input}/*{1,2}.fastq"


// PRINT STANDARD LOGGING INFO
include { printLogging } from './lib/functions.nf'
printLogging(reads_path)



////////////////////
// STAGE CHANNELS //
////////////////////

/*
 *   Channels are where you define the input for the different
 *    processes which make up a pipeline. Channels indicate
 *    the flow of data, i.e. the "route" that a file will take.
 */

// STAGE BAM FILES FROM TEST PROFILE # this establishes the test data to use with -profile test
if ( workflow.profile.tokenize(',').contains('test') ){

        include { check_test_data } from './lib/functions.nf'
        READS = check_test_data(params.readPaths, params.SE)

} else {

    // STAGE READS CHANNELS # this defines the normal input when test profile is not in use
    READS = Channel
        .fromFilePairs(reads_path, size: params.SE ? 1 : 2)
        .ifEmpty{ exit 1, """Cannot find valid read files in dir: ${params.input}
        The pipeline will expect PE reads in compressed *{1,2}.${params.extension} format
        unless you have specified the --SE parameter or a different extension using --extension """}
        .map{ tuple(it[0], it[1]) }
        .take(params.take.toInteger())

}



////////////////////
// BEGIN WORKFLOW //
////////////////////

/*
 *   Workflows are where you define how different processes link together. They
 *    may be modularised into "sub-workflows" which must be named eg. 'DNAseq'
 *    and there must always be one MAIN workflow to link them together, which
 *    is always unnamed.
 */

include { DNAseq } from "${projectDir}/modules/workflow"

// MAIN WORKFLOW 
workflow {

    // call sub-workflows eg. WORKFLOW(Channel1, Channel2, Channel3, etc.)
    main:
        DNAseq(READS, fasta, fai)

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
include { printSummary } from './lib/functions.nf'
workflow.onComplete {

    printSummary()

    // run a small clean-up script to remove "work" directory after successful completion 
    if (!params.debug && workflow.success) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}
