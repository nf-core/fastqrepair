<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-fastqrepair_logo_dark.png">
    <img alt="nf-core/fastqrepair" src="docs/images/nf-core-fastqrepair_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/fastqrepair/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/fastqrepair/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/fastqrepair/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/fastqrepair/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/fastqrepair/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

<!-- [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) -->

[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/fastqrepair)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23fastqrepair-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/fastqrepair)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/fastqrepair** is a bioinformatics pipeline that can be used to recover corrupted `FASTQ.gz` files, drop or fix uncompliant reads, remove unpaired reads, and settles reads that became disordered. It takes a `samplesheet` and FASTQ/FASTQ.gz files as input (both single-end and paired-end) and produces clean FASTQ files and a QC report.

![pipeline_diagram](docs/images/fastqrepair-flow-diagram-v1.0.svg)

1. Recover reads from corrupted fastq.gz file ([`gzrt`](https://github.com/arenn/gzrt))
2. Make recovered reads well-formed ([`fastqwiper`](https://github.com/mazzalab/fastqwiper))
3. Drop unpaired reads ([`trimmomatic`](http://www.usadellab.org/cms/index.php?page=trimmomatic))
4. Re-pair reads ([`bbmap/repair.sh`](https://sourceforge.net/projects/bbmap/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

**samplesheet.csv**:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
mysampleA,sample_R1.fastq.gz,sample_R2.fastq.gz
mysampleB,sample_R3.fastq.gz,sample_R4.fastq.gz
mysampleC,sample_R5.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end). Rows with the same sample identifier are not allowed. Row with different sample identifiers but same file names are not allowed.

Now, you can run the pipeline using:

```bash
nextflow run nf-core/fastqrepair \
   -profile <test/docker> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

optional parameters are:

```bash
--chunk_size  <int multiple of 4>
--qin         <33/64>
--alphabet    <ACGTN>
```

where

`chunk_size` is the number of lines of chunks of the original fastq file (caution! Too big or too small numbers may significantly impact on performance); `qin` is the ASCII offset (33=Sanger, 64=old Solexa); `alphabet` is the allowed alphabet in the SEQ line of the FASTQ file.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/fastqrepair/usage) and the [parameter documentation](https://nf-co.re/fastqrepair/parameters).

## Pipeline output

This pipeline produces clean and well-formed fastq files together with short textual reports of the cleaning actions.

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/fastqrepair/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/fastqrepair/output).

## Credits

`nf-core/fastqrepair` was designed and written by [Tommaso Mazza](https://github.com/mazzalab).

<!-- We thank the following people for their extensive assistance in the development of this pipeline: -->
<!-- nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#fastqrepair` channel](https://nfcore.slack.com/channels/fastqrepair) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/fastqrepair for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
