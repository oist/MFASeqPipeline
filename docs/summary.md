# Workflow engine summary

`arubaito.py` is a small, self-contained workflow engine (~600 lines, only
dependency PyYAML). A run is driven by one YAML *configuration* which defines
**variables**, **steps** and **targets**, and may pull in other YAML files via
**`use`**. Each step's shell script is launched through a configurable
**executor** (`execute_cmd`), e.g. directly or via Slurm `srun` on an HPC
cluster, with optional staging of input/output files through shared temporary
storage. The engine has no knowledge of MFA-seq; all domain logic lives in the
YAML files and the helper scripts under `steps/`.

```
arubaito [options] CONFIG [TARGET ...] [KEY=VALUE ...]      # installed with pip
./arubaito.py [options] CONFIG [TARGET ...] [KEY=VALUE ...] # from a checkout
```

## 1. Configuration files and `use`

* `CONFIG` is a name (`example_dataset`) or file (`example_dataset.yaml`). Files
  are searched in `./` (the current directory) first, then in the directory
  containing `arubaito.py` (for a pip install, the installed package
  `mfaseq_pipeline`, which also contains all pipelines; see section 8). The same
  search applies to every `use` entry.
* `use: name` or `use: [name1, name2]` loads other configs recursively and merges
  them with a **shallow, top-level** dict update:
  * The including file always wins over what it uses.
  * In a list, **earlier entries win over later ones** (`use: [a, b]`: a value
    defined in both `a` and `b` comes from `a`).
  * Shallow means a top-level key is replaced wholesale — e.g. redefining
    `map_step` replaces the whole step, redefining `targets` replaces all targets.
* For every loaded file `<name>.yaml` the engine sets `<name>_dir` to the search
  directory where it was found (e.g. `mfa_pipeline_dir`). Pipelines use this to
  call their helper scripts: `python3 {mfa_pipeline_dir}/steps/abundance.py`. Note that
  for files found in the current directory this is the relative path `./`.

The engine does not distinguish kinds of YAML files: any file may define
variables, steps and targets, and any file may `use` any other. The usual
*convention* (see `example_dataset.yaml`) is a *project/dataset* file that sets
paths, parameters and `targets`, and `use`s one or more reusable *pipeline* files
(`*_pipeline.yaml`) that define default variables and steps. But a workflow may
just as well define its own steps directly, with no separate pipeline file, and
a pipeline file could define targets.

## 2. Variables

Every top-level key that isn't otherwise special is a variable. Values are
usually strings and may reference other variables as `{name}`.

* Expansion is Python `str.format_map`, applied **repeatedly until the string
  stops changing** (max ~100 rounds; this detects cycles). So variables can refer
  to variables that refer to variables, etc.
* `{{` and `}}` yield literal `{` and `}` (the engine protects them across the
  repeated rounds). This matters mostly in shell scripts: `${{var}}`,
  `"${{arr[@]}}"`, awk `{{ print $2 }}`.
* Expansion is **lazy / late-bound**: values are stored unexpanded and expanded
  only when used, in the context of the current job and step. Thus a global
  `bam_file: '{out_dir}/{sample}.bam'` works even though `sample` is only
  defined per target.
* Undefined references raise a `ConfigurationError` listing the expansion chain.
* Non-string values (numbers, booleans, lists) are inserted via `str()`.
* The engine does pure text substitution. Things like
  `step_memory_gb: '$[{memory_map_gb}+{memory_sort_gb}]'` or
  `growthrate_perhour: $(awk ...)` only work because the result ends up inside a
  bash script or `execute_cmd` and is evaluated by the shell.

Only these things are ever expanded: step scripts, input/output file variables,
`log_file`, `execute_cmd`, and `tmp_dir`.

### Lookup precedence (highest first)

| Layer | Source |
|---|---|
| step-internal | `step` (current step name), resolved input/output paths, `script_file` |
| step | the step's `override:` dict |
| command line | `KEY=VALUE` arguments and `--steps` |
| job | the job's own variables in its target list |
| target groups | `override:` dicts of enclosing groups, innermost first |
| global | merged top-level config |

Note that a step's `override` beats even the command line (e.g. `mapstats` forces
`step_cores: 1`).

Input/output file variables are expanded in the **job** context, i.e. *without*
the step override. This is deliberate: files are how steps communicate, so the
output path one step writes must be exactly the input path the next step reads,
whatever per-step settings each of them has. Scripts, `log_file` and
`execute_cmd` are expanded in the **step** context. (The path a script actually
sees may still differ from that canonical path when temp staging is on — see
below.)

## 3. Steps

A step is a top-level key `<name>_step`:

```yaml
abundance_step:
  input:  [bam_file, mappability_file, reference_bam_file]   # variable NAMES
  output: [abundance_file]                                   # variable NAMES
  override:                                                  # optional
    step_cores: 1
    log_label: '.w{windowsize}bp'
  script: |
    python3 {mfa_pipeline_dir}/steps/abundance.py --sample-bam "{bam_file}" ... --output "{abundance_file}"
```

* `input` / `output` list **variable names**, not paths. Each is expanded to its
  canonical (final) path. Inside the script, `{bam_file}` is rebound to the path
  the script should actually use: the canonical path, or a temporary copy when
  temp staging is enabled (see below). Scripts must therefore always refer to
  their files through these variables, never by reconstructing paths.
* An input/output that expands to `''` is ignored — this is how optional files
  are switched off (e.g. `mfa_no_reference.yaml` blanks the reference files;
  `velogp_model_file: ''`).
* **Skip rule:** a step is skipped iff *all* its (non-empty) outputs exist.
  There is no timestamp/staleness check; delete outputs to force a rerun.
  Consequence: **a step with no outputs is always skipped**.
* Inputs are not checked for existence, and there is no dependency resolution.
  Steps run strictly in the order of the `steps` list; `input`/`output` only
  serve the skip rule and temp-dir staging.
* After the script runs, missing outputs cause an `ExecutionError`.

### How a step is executed

The engine never runs step scripts directly; it hands them to an **executor**,
a command template in the variable `execute_cmd`. This is what lets the same
pipeline run locally or on an HPC cluster.

1. The script is expanded (step context) and written to a temp file in the job's
   `tmp_dir`, with shebang `#!{shell}`; `script_file` is set to its path.
2. `execute_cmd` is expanded (step context) and run through `{shell}`; the
   engine waits for it to finish and treats a non-zero exit as failure.
   * Locally: `execute_cmd: '{script_file}'`.
   * On a Slurm cluster (default in `mfa_pipeline.yaml`):
     `srun -p '{step_partition}' --job-name=... --mem={step_memory_gb}G --time={step_time} -c {step_cores} {script_file}`.
     The per-step resource variables (`step_partition`, `step_cores`,
     `step_memory_gb`, `step_time`) exist for this purpose and are tuned via step `override`s.
   * Because the script then runs on a compute node, `tmp_dir` (which holds the
     script file) must be visible from there.
3. stdout/stderr go to `log_file` if defined (expanded per step; default in
   `mfa_pipeline.yaml`: `{out_dir}/{target_label}.{step}{log_label}.log`), and
   are additionally `tee`d to the terminal unless `--quiet`.
4. Each log starts with a header listing the step, its inputs, outputs and log.

### Temp-dir staging (cluster support)

On a cluster, the storage holding inputs and results is not necessarily
accessible from the compute nodes, or only read-only. The engine therefore
distinguishes a file's canonical path from the path a step script uses, and can
stage files through the job's `tmp_dir`, which must be on storage visible to
*both* the machine running the engine and the compute nodes (shared-nothing
clusters are not supported):

* Each job gets a private directory from `mkdtemp`, using `tmp_dir` as *prefix*
  (so `/tmp/` → `/tmp/<random>`, `/tmp/x-` → `/tmp/x-<random>`). `tmp_dir` is then
  rebound to the actual directory and deleted when the job ends.
* `execute_copy_inputs: True` — inputs are copied into `tmp_dir` (files or whole
  directories) and scripts see the copies. Inputs from the same source directory
  share a temp subdirectory (keeps e.g. a BAM and its index together).
* `execute_copy_outputs: True` — scripts write outputs into `tmp_dir`; the
  engine copies them to their final location afterwards. Later steps of the same
  job read the temp copy. (Default in `mfa_pipeline.yaml`, for result storage
  that is read-only on compute nodes.)
* The skip rule and the existence check after a step always use the canonical
  paths.
* Both flags may also be set per step via `override`.

## 4. Targets and jobs

The top-level `targets` key defines what to run. The unit of execution is a
**job**: a dict of variable bindings. A **target** is a named list of jobs.

```yaml
targets:
  stationary:                      # group as a plain list of jobs
    - sample: s1
    - sample: s2
  exponential:                     # group as a dict
    override:                      # variables shared by all jobs below
      steps:
        - map
        - mapstats
        - abundance
    targets:                       # the group's own job list
      - {sample: e1, reference: s1}
      - {sample: e2, reference: s2}
    late:                          # any other key = named sub-target
      - {sample: e3, reference: s1}
```

* A target value is either a **list** (each element is one job dict) or a
  **dict** with optional special keys `override` (shared variables) and `targets`
  (job list or nested group); all other keys are named sub-targets.
* Target names are a **single flat namespace**: every name at every nesting level
  becomes a target that contains all jobs below it. Above, `exponential` has 3
  jobs and `late` 1. Equal names in different places get merged.
* The special variable **`steps`** selects which steps a job runs (default in
  `mfa_pipeline.yaml`). It is an ordinary variable, so it can be set globally, per
  group, per job, or via `--steps a,b,c`.
* Convention: `target_label` (default `{sample}`) names logs, tmp dirs and
  cluster jobs.
* Elements of a job list are always taken as job dicts; a list element of the
  form `{override: ..., targets: ...}` is *not* treated as a nested group. Use a
  dict with named sub-targets instead (see the third example in the comment
  above `add_target`).

### Running targets

* Targets named on the command line run **sequentially in the given order**;
  jobs within one target run in parallel with `--jobs N` (thread pool).
  So dependencies between jobs (e.g. reference samples, genome index) are handled
  by putting them in separate targets and listing those first.
* `all` expands to all **top-level** target names.
* With no target, a single job runs with only the global configuration.
* Everything from the first `KEY=VALUE` argument on is a setting override;
  values are parsed as YAML (`windowsize=5000`, `steps=[map]`).

### Options

Options must come **before `CONFIG`** or **after everything else**. An option
placed between `CONFIG`, the targets and the `KEY=VALUE` settings splits that
list, and the arguments after it are rejected as "unrecognized arguments"
(e.g. `arubaito.py cfg t1 --dry t2` fails; `arubaito.py --dry cfg t1 t2` works).

| Option | Effect |
|---|---|
| `--list-targets` | print every target with its jobs' variables, exit |
| `--jobs N` | parallel jobs per target (implies `--quiet`) |
| `--quiet` / `--no-quiet` | don't echo logs to stdout |
| `--dry` | print step headers, do not run scripts (still creates tmp/output dirs) |
| `--catch` / `--no-catch` | continue with other jobs after a failure (default on when targets are given) |
| `--keep` | on failure, still copy produced outputs out of `tmp_dir` |
| `--verbose` | log copies, the expanded script and the executor command |
| `--steps a,b` | override `steps` |

## 5. Engine-level variables

| Variable | Role |
|---|---|
| `steps` | list of step names to run for a job |
| `targets` | target definitions |
| `use` | configs to include |
| `shell` | interpreter for scripts and `execute_cmd` (`/bin/bash`) |
| `execute_cmd` | executor: how to launch `{script_file}` (locally, `srun`, …) |
| `tmp_dir` | prefix for the per-job temp directory (rebound to the real dir) |
| `execute_copy_inputs`, `execute_copy_outputs` | temp-dir staging flags |
| `log_file` | optional per-step log path |
| `step`, `script_file`, `<name>_dir` | set automatically |

Everything else (`step_partition`, `step_cores`, `step_memory_gb`, `step_time`, `target_label`,
`log_label`, `sample`, `out_dir`, …) is pure convention of the pipeline YAMLs.

## 6. Bundled configuration files

| File | Content |
|---|---|
| `mfa_pipeline.yaml` | main paired-end MFA-seq pipeline: variables and steps `genome_index`, `map`, `fastqc`, `mapstats`, `mappability`, `abundance`, `abundance_plot`, `velocity_oscillation`, `velogp`, `velomix` |
| `mfa_se_pipeline.yaml` | single-end variant; `use: mfa_pipeline`, replaces `map`/`mappability`, `read_file`, `paired_end: False` |
| `mfa_no_reference.yaml` | overlay that blanks the `reference_*` inputs |
| `rnaseq_pipeline.yaml` | `genecounts` step (featureCounts); no `use`, relies on `mfa_pipeline` variables |
| `refseq_pipeline.yaml` | NCBI genome/GFF/GBFF download steps (`genome_download`, `genome_gbff_download`, `genome_download_faonly`) |
| `sra_pipeline.yaml`, `sra_se_pipeline.yaml` | SRA download via `prefetch` + `fastq-dump` (`sra`, `fastq` steps) |
| `sra_pipeline2.yaml`, `sra_se_pipeline2.yaml` | `fastq` step using `fasterq-dump` directly; `accession` may be comma-separated |
| `example_dataset.yaml` | template for a project config with targets |
| `env.yaml` | conda environment |

Helper programs live in `steps/` (`abundance.py`, `mappability.py`,
`mapping_statistics.py`, R scripts) and the git submodule `steps/velogp`.

## 7. Known quirks

* A step without outputs never runs (skip rule on an empty list).
* No staleness detection: changed inputs or parameters do not trigger reruns
  unless the parameter is part of an output file name (e.g. `abundance_file`
  contains `{windowsize}`).
* With `execute_copy_outputs`, two steps in one job producing the same output
  variable raise an error.

## 8. Packaging (pip)

`pyproject.toml` makes the repository root the Python package `mfaseq_pipeline`
(distribution name `mfaseq-pipeline`, `package-dir = {"mfaseq_pipeline" = "."}`,
plus an empty root `__init__.py`). The engine, all `*.yaml` files (including
`env.yaml`) and everything under `steps/` (including the VeloGP submodule) are
installed together into `site-packages/mfaseq_pipeline/`, so the engine finds the
pipelines next to itself and `{<name>_dir}/steps/...` resolves exactly as in a
checkout. The console script `arubaito` calls `mfaseq_pipeline.arubaito:main`.

* Helper scripts are invoked through their interpreter (`python3 …`, `Rscript …`)
  rather than relying on the executable bit.
* pip cannot install the external tools; users create the conda environment from
  `env.yaml` (downloaded from GitHub, or from a clone) and pip-install into it.
* `pip install -e .` from a clone uses the clone's files directly.
* Building creates `mfaseq_pipeline.egg-info/` in the repository root; it is
  covered by `.gitignore`.
