# nf-core/fastqrepair: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v[1.1.1](https://github.com/nf-core/fastqrepair/releases/tag/1.1.1) - Reggina Amaranth [27/08/2026]

First release since v1.0.0. Includes the never-released v1.1.0 changes plus a nf-core/tools v4.0.3 template update (supersedes the pending v3.4.1, v3.5.1 and v4.0.2 template merges).

### `Changed`

- [PR #22](https://github.com/nf-core/fastqrepair/pull/22) - Template update for nf-core/tools v3.3.2
- [#17](https://github.com/nf-core/fastqrepair/issues/17) - Improved `nextflow_schema.json` setting a minimum of 1 split to `num_splits` and removing a regex
- Template update for nf-core/tools v4.0.3. Minimum Nextflow version raised to `25.10.4`; `nf-schema` bumped to `2.5.1`; `MultiQC` bumped to `1.34`
- `FASTQREPAIR` now takes `multiqc_config`, `multiqc_logo`, `multiqc_methods_description` and `outdir` as explicit workflow inputs, and collates software versions through the `versions` channel topic

### `Fixed`

- `FASTQC.out.versions` no longer exists after the template's `FASTQC` module switched to publishing its version through the `versions` topic channel instead of a plain `emit: versions` output; the stale reference in `workflows/fastqrepair.nf` crashed every run with a `MissingPropertyException`
- `tests/main.nf.test` expected the pre-4.0 collated-versions filename (`nf-core_fastqrepair_versions.yml`); updated to the template's new name (`nf_core_fastqrepair_software_mqc_versions.yml`)
- `tests/.nftignore` patterns assumed the nf-core default `multiqc/`/`fastqc/` output layout; this pipeline publishes under `QC/multiqc/`/`QC/fastqc/`, so the patterns never matched and non-deterministic files (FastQC's zip, MultiQC's HTML report and SVG plots) were being hashed into `default.nf.test`'s snapshot, making it fail on every rerun

### `Removed`

- `--hook_url` parameter and the Slack/Microsoft Teams notification assets (`assets/slackreport.json`, `assets/adaptivecard.json`), dropped by the template
- `gitpod` profile and `.gitpod.yml`; the `arm` profile is replaced by `arm64` plus `emulate_amd64`

## v[1.0.0](https://github.com/nf-core/fastqrepair/releases/tag/1.0.0) - Catanzaro YellowRed [04/02/2025]

Initial release of nf-core/fastqrepair, created with the [nf-core](https://nf-co.re/) template.

### `Fixed`

- [PR #2](https://github.com/nf-core/fastqrepair/pull/2) - First release

### `Dependencies`

| Dependency   | Old version | New version |
| ------------ | ----------- | ----------- |
| `BBmap`      |             | 39.13       |
| `FastQC`     |             | 0.12.1      |
| `gzrt`       |             | 0.9.1       |
| `Wipertools` |             | 1.1.5       |
| `MultiQC`    |             | 1.26        |

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

### Credits

Special thanks to the following for their contributions to the release:

- [Sateesh Peri](https://github.com/sateeshperi)
- [Louis Le Nézet](https://github.com/LouisLeNezet) - reviewer
- [Anabella Trigila](https://github.com/atrigila) - reviewer
- [James A. Fellows Yates](https://github.com/jfy133) - reviewer
- [Charles Plessy](https://github.com/charles-plessy) - reviewer
