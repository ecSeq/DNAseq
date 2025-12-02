// INCLUDES # here you must give the relevant processes from the modules/process directory 
include {cutadapt;FastQC;bowtie2_index;bowtie2;bamQC} from "${projectDir}/modules/process"

// SUB-WORKFLOWS
workflow 'DNAseq' {

    // take the initial Channels and paths
    take:
        READS
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
        bowtie2_index(fasta,fai)

        // now we have a series of trimmed reads from cutadapt eg. [sample_name, [clipped_1.fastq.gz,clipped_2.fastq.gz]]
        // we also have a series of indexed genomes for mapping eg. [/path/to/index/dir]
        // so we should combine these channels to make all possible combinations:
        reads_and_index = cutadapt.out[0].combine(bowtie2_index.out)
        // eg. [sample_name, [clipped_1.fastq.gz,clipped_2.fastq.gz], /path/to/index/dir]

        // and now we can run the alignments!
        bowtie2(reads_and_index,fasta,fai)

        // finally we are able to run the optional bamQC if the user chooses to
        bamQC(bowtie2.out[0])

}
