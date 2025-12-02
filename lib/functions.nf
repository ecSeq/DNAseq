// This file for loading custom functions into the main.nf script (separated for portability)

// print --help information
def printHelp() {
    println """\

     ===========================================
      E C S E Q - D N A s e q   P I P E L I N E
     ===========================================
     ~ version ${workflow.manifest.version}

     Usage: 
          nextflow run ecseq/dnaseq [OPTIONS]...

     Options: GENERAL
          --input [path/to/input/dir]     [REQUIRED] Provide the directory containing
                                           fastq file(s) in "*{1,2}.fastq.gz" format

          --reference [path/to/ref.fa]    [REQUIRED] Provide the path to the reference
                                           genome in fasta format

          --output [STR]                  A string that can be given to name the output
                                           directory. [default: "."]

          --SE                            Indicate to the pipeline whether fastq files
                                           are SE reads in "*.fastq.gz" format.
                                           [default: off]


     Options: MODIFIERS
          --FastQC                        Generate FastQC report of trimmed reads.
                                           [default: off]

          --bamQC                         Generate bamQC report of alignments.
                                           [default: off]

          --keepReads                     Keep trimmed fastq reads. [default: off]


     Options: TRIMMING
          --forward                       Forward adapter sequence.
                                           [default: "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"]

          --reverse                       Reverse adapter sequence.
                                           [default: "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"]

          --minQual                       Minimum base quality threshold. [default: 20]

          --minLeng                       Minimum read length threshold. [default: 25]

          --minOver                       Minimum overlap threshold. [default: 3]


     Options: ADDITIONAL
          --help                          Display this help information and exit
          --version                       Display the current pipeline version and exit
          --debug                         Run the pipeline in debug mode    


     Example: 
          nextflow run ecseq/dnaseq \\
           --input /path/to/input/dir \\
           --reference /path/to/genome.fa \\
           --FastQC --bamQC --keepReads

    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// print --version information
def printVersion() {
    println """\
         ===========================================
          E C S E Q - D N A s e q   P I P E L I N E
         ===========================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// print pipeline initiation logging
def printLogging(reads_path) {
    log.info """
         ===========================================
          E C S E Q - D N A s e q   P I P E L I N E """
    if(params.debug){ log.info "         (debug mode enabled)" }
    log.info """         ===========================================
         ~ version ${workflow.manifest.version}

         input dir    : ${workflow.profile.tokenize(',').contains('test') ? '-' : reads_path}
         reference    : ${params.reference}
         output dir   : ${params.output}
         mode         : ${params.SE ? 'single-end' : 'paired-end'}
         QC options   : ${params.FastQC ? 'FastQC ' : ''}${params.bamQC ? 'bamQC' : ''}

         ===========================================
         RUN NAME: ${workflow.runName}

    """
}

// print pipeline execution summary
def printSummary() {

    log.info """
         Pipeline execution summary
         ---------------------------
         Name         : ${workflow.runName}${workflow.resume ? ' (resumed)' : ''}
         Profile      : ${workflow.profile}
         Launch dir   : ${workflow.launchDir}
         Work dir     : ${workflow.workDir} ${params.debug || !workflow.success ? '' : '(cleared)' }
         Status       : ${workflow.success ? 'success' : 'failed'}
         Error report : ${workflow.errorReport ?: '-'}
    """
}

// FUNCTION TO LOAD DATASETS IN TEST PROFILE
def check_test_data(readPaths, singleEnd) {

    // Set READS testdata
    if( singleEnd ) {
        READS = Channel.from(readPaths)
            .transpose()
            .map { row -> [ file(row[1]).getName().tokenize(".").get(0), file(row[1]) ] }
            .ifEmpty { exit 1, "test profile readPaths was empty - no input files supplied" }
    } else {
        READS = Channel.from(readPaths).view()
            .map { row -> [ row[0], [file(row[1][0]),file(row[1][1])] ] }
            .ifEmpty { exit 1, "test profile readPaths was empty - no input files supplied" }
    }
    // Return READS channel
    return READS
}
