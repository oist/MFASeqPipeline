# refseq_pipeline

File: `refseq_pipeline.yaml`.

## Summary

`refseq_pipeline` downloads a reference genome and its annotation from NCBI
using the [NCBI Datasets API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/),
given an assembly accession (variable `refseq`, e.g. `GCF_000005845.2` for
*E. coli* K-12 MG1655). The files are stored as `{refseq_dir}/{refseq}.fna` /
`.gff` / `.gbff`, with `refseq_dir` defaulting to `{project_dir}/refseq`. Because
it defines `genome_file` and `gff_file`, it can be used together with
[mfa_pipeline](mfa_pipeline.md) or [rnaseq_pipeline](rnaseq_pipeline.md) to
analyse data against an NCBI genome. Downloads need `curl`, `unzip` and network
access from wherever the steps run.

## Quick start

Stand-alone, the pipeline only needs a few engine settings. Create `genomes.yaml`:

```yaml
use: refseq_pipeline
project_dir: /data/myproject
tmp_dir: /tmp/refseq-
shell: /bin/bash
execute_cmd: '{script_file}'
targets:
  ecoli:
    - refseq: GCF_000005845.2
      steps: [ genome_download, genome_gbff_download ]
```

```bash
./arubaito.py genomes ecoli
```

This creates `/data/myproject/refseq/GCF_000005845.2.fna`, `.gff` and `.gbff`.

Combined with `mfa_pipeline`, list `refseq_pipeline` **first** in `use`. Earlier
entries take precedence, so its `genome_file` and `gff_file` replace the
defaults of `mfa_pipeline`:

```yaml
use: [ refseq_pipeline, mfa_pipeline ]
refseq: GCF_000005845.2
targets:
  genome:
    - target_label: genome
      steps: [ genome_download, genome_index ]
  # ... sample targets as in the mfa_pipeline quick start
```

## Steps

All steps use one core and 1 GB of memory, and set `target_label` to `{refseq}`.

### genome_download

Downloads the genome sequence and annotation.

* Setting `refseq`: NCBI assembly accession (`GCF_…` or `GCA_…`).
* **Output** `genome_file` (`{refseq_dir}/{refseq}.fna`): genome sequence,
  [FASTA](https://en.wikipedia.org/wiki/FASTA_format). Sequence names are NCBI
  accessions such as `NC_000913.3`.
* **Output** `gff_file` (`{refseq_dir}/{refseq}.gff`): genome annotation,
  [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).

### genome_download_faonly

Same as `genome_download`, but keeps only the sequence. Use this for assemblies
without an annotation.

* **Output** `genome_file` (`{refseq_dir}/{refseq}.fna`): FASTA.

### genome_gbff_download

Downloads the annotated genome in GenBank format.

* **Output** `gbff_file` (`{refseq_dir}/{refseq}.gbff`):
  [GenBank flat file](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
  (sequence plus annotation in one file).
