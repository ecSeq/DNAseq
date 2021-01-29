[<img width="200" align="right" src="docs/images/ecseq.jpg">](https://www.ecseq.com)
[![Nextflow](https://img.shields.io/badge/nextflow-20.07.1-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/ecseq/dnaseq.svg)](https://hub.docker.com/r/ecseq/dnaseq)

ecSeq-DNAseq Pipeline
======================

**ecSeq/DNAseq** is a simple bioinformatics analysis pipeline for trimming and aligning DNAseq data with a small selection of software.

The workflow processes a collection of raw fastq files using [cutadapt](https://github.com/marcelm/cutadapt), producing quality trimmed reads which are then evaluated with [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Indexing of the given reference genome is performed prior to alignment with any combination of [bowtie2](https://github.com/BenLangmead/bowtie2), [BWA/BWA MEM](https://github.com/lh3/bwa), [segemehl](https://www.bioinf.uni-leipzig.de/Software/segemehl/) and [STAR](https://github.com/alexdobin/STAR). 

> See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://www.nextflow.io/)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run ecseq/dnaseq -profile test,<docker|singularity|conda>
```

iv. Start running your own analysis!

```bash
nextflow run ecseq/dnaseq -profile <docker|singularity|conda> --input /path/to/fastq/dir --reference /path/to/genome.fa
```

> See the [usage documentation](docs/usage.md) for all of the available options when running the pipeline.


### Credits

These scripts were originally written for use by [ecSeq Bioinformatics GmbH](https://www.ecseq.com), by Adam Nunn ([@bio15anu](https://github.com/bio15anu)).
