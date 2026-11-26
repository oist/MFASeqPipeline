# arubaito

`arubaito` is the small workflow engine that runs the MFASeqPipeline pipelines.
It is a single Python script, `arubaito.py`, configured by YAML files. Installing
MFASeqPipeline with pip (see the [README](../README.md#installation)) provides it
as the command `arubaito`. Without installing, run `arubaito.py` from a checkout
of the repository instead; everything described here works the same way.

arubaito executes **jobs**. A job is a sequence of **steps**: map the reads,
compute statistics, compute abundances, and so on. Each job has its own set of
**variables**, e.g. the sample it processes. Jobs are grouped into
**targets**, and you run a workflow by naming the targets to execute:

```bash
arubaito myproject stationary exponential
```

This loads `myproject.yaml`, then runs all jobs of target `stationary`, then all
jobs of target `exponential`.

## Steps

A step is a shell script with placeholders. Before the script runs, arubaito
substitutes the job's variables into it: `{sample}` becomes the sample name,
`{bam_file}` the path of the BAM file, and so on. Substitution is **recursive**:
if the value of a variable contains `{...}` placeholders itself, those are
substituted too, until none are left.

A step named `map` is defined by the key `map_step`:

```yaml
map_step:
  input:  [ read1_file, read2_file, genome_index_dir ]
  output: [ bam_file, bam_index_file ]
  override:
    step_memory_gb: 48
  script: |
    bwa mem -t {step_cores} "{genome_index_prefix}" "{read1_file}" "{read2_file}" \
      | samtools sort -o "{bam_file}" -
    samtools index "{bam_file}" "{bam_index_file}"
```

* `input` and `output` list **variable names**, not file names. Each of them names a
  variable holding a file path.
* `override` (optional) sets variables for this step only.
* `script` is the shell script (bash by default).

A step is **skipped if all its outputs already exist**. To redo a step, delete its
outputs. There is no check whether inputs are newer than outputs. Steps that have
no outputs at all are therefore always skipped. If a step finishes without
creating all its outputs, the job fails. An output variable whose value is empty
(`''`) is ignored, which is how optional outputs are switched off.

The variable `steps` lists the steps a job runs, in that order. Like any other
variable, it can be set globally, per target or per job, or on the command line
with `--steps map,mapstats`. Steps run in exactly the listed order; arubaito does
not work out dependencies between steps.

## Variables

Variables are the top-level keys of the configuration (other than the step
definitions, `targets` and `use`). Most values are strings that may contain
placeholders:

```yaml
project_dir: '{base_dir}/{project}'
out_dir: '{project_dir}/results'
bam_file: '{out_dir}/{sample}.bam'
```

### Expansion rules

* `{name}` is replaced by the value of variable `name`, and the result is expanded
  again, until nothing changes. Above, `{bam_file}` becomes
  `/data/myproject/results/exp1.bam`.
* `{{` and `}}` stand for literal `{` and `}`. This is needed for shell syntax in
  scripts, e.g. `"${{files[@]}}"` or `awk '{{ print $2 }}'`.
* Expansion happens **when a step runs**, using the variables of that job and step.
  So a global setting like `bam_file` may refer to `{sample}`, which only targets
  define.
* Referring to an undefined variable is an error, and so is a circular
  definition (`a: '{b}'`, `b: '{a}'`).
* arubaito only substitutes text. Shell syntax in values, like
  `$(awk ...)` or `$[1 + 2]`, is evaluated by bash when the script runs.
* Numbers and booleans are inserted as text (`True`, `4`).

### Where variables come from

When several places define the same variable, the first one in this list wins:

1. the step's `override`
2. `KEY=VALUE` settings on the command line, and `--steps`
3. the job's own variables
4. the `override` of the target groups enclosing the job, innermost first
5. the configuration file and the files it `use`s

A step's `override` takes precedence even over the command line. For example,
the `mapstats` step always uses one core, whatever `cores` is set to.

There is one exception: **input and output paths** are expanded *without* the
step's `override`. Files are how steps pass data to each other, so the file one
step writes must have exactly the name the next step reads, whatever settings
either step overrides.

### Configuration files and `use`

A configuration can include other files with `use`:

```yaml
use:
  - refseq_pipeline
  - mfa_pipeline
```

* Files are looked up in the current directory first, then next to `arubaito.py`.
  When MFASeqPipeline is installed with pip, that is the installed package,
  which contains all pipelines.
* The including file wins over the files it uses, and **earlier entries in the
  list win over later ones**. Above, `refseq_pipeline`'s `genome_file` replaces
  `mfa_pipeline`'s.
* Merging happens at the top level only: redefining `map_step` replaces the
  whole step, and redefining `targets` replaces all targets.
* For each loaded file `<name>.yaml`, the variable `<name>_dir` holds the directory
  it was found in, e.g. `mfa_pipeline_dir`.

## Targets

Targets are defined under the top-level key `targets`. In the simplest case a
target is a list of jobs, and each job is a set of variables:

```yaml
targets:
  stationary:
    - sample: stat1
    - sample: stat2
```

### Nesting and overrides

A target can also be a **group** (a dictionary) that contains:

* `override`: variables shared by all jobs in the group;
* `targets`: the group's own list of jobs;
* any other key: a named sub-target, which can again be a list or a group.

```yaml
targets:
  experiment:
    override:                       # applies to all jobs below
      steps:
        - map
        - mapstats
        - mappability
        - abundance
        - velogp
      windowsize: 5k
    stationary:
      - sample: stat1
      - sample: stat2
    exponential:
      override:                     # applies to exp1 and exp2 only
        windowsize: 10k
      targets:
        - { sample: exp1, reference: stat1 }
        - { sample: exp2, reference: stat2, windowsize: 20k }
```

* Every named group, at any depth, is a target of its own containing all jobs
  below it. Here `experiment` has four jobs, `stationary` two and `exponential` two.
  Target names are global: a name used in two places refers to the jobs of both.
* Variables are looked up from the inside out: the job first, then the innermost
  group's `override`, then the next group out. So `exp2` uses `windowsize: 20k`,
  `exp1` uses `10k`, and `stat1` and `stat2` use `5k`.
* The entries of a plain list are always jobs. To nest groups, use a dictionary
  with named sub-targets, as above.

### Running targets

* Targets named on the command line run **one after another**, in the order given.
  Put targets that others depend on first, e.g. the reference samples before the
  samples that use them.
* Jobs **within** a target are independent and can run in parallel:
  `--jobs 4` runs up to four at a time.
* `all` stands for all top-level targets.
* `KEY=VALUE` arguments after the targets override variables for this run:
  `arubaito myproject exponential windowsize=20k`.
* Options must come before the configuration name or at the very end of the
  command line. `--list-targets` shows all targets and their jobs. `--dry` shows
  what would run without running it.

## HPC clusters

### The executor

arubaito never runs a step's script directly. It writes the expanded script to a
file and runs the command in the variable **`execute_cmd`**, with `{script_file}`
replaced by that file's path. To run steps locally:

```yaml
execute_cmd: '{script_file}'
```

To submit each step to a Slurm cluster (the default in `mfa_pipeline`):

```yaml
execute_cmd: srun -p '{step_partition}' --job-name='{target_label}.{step}' --mem={step_memory_gb}G --time={step_time} -c {step_cores} {script_file}
```

* The command must **wait** until the step has finished and return its exit
  status, as `srun` does. arubaito then checks the outputs and starts the next
  step.
* `execute_cmd` is expanded like a script, so it can use per-step variables. In
  the pipelines, `step_partition`, `step_cores`, `step_memory_gb` and `step_time`
  exist for this purpose. Steps set them in their `override`; for example, `map`
  asks for more memory.
* Use `--jobs N` to have up to N cluster jobs running at once.
* The step's output (stdout and stderr) is written to the file in `log_file`, by
  default `{out_dir}/{target_label}.{step}.log`.
* The software the steps need, e.g. the conda environment, must be available on
  the compute nodes. `srun` passes the current environment on by default.

### Temporary storage and copying

Each job gets its own temporary directory, created with the value of `tmp_dir`
as a name prefix and deleted when the job ends. The script files are placed
there, so **`tmp_dir` must be on storage that the compute nodes can access**.

On some clusters, the storage holding the input data or the results is not
available on the compute nodes, or only read-only. Two settings let steps work
on copies in `tmp_dir` instead:

* `execute_copy_outputs: True`: steps write their outputs into `tmp_dir`, and
  arubaito copies them to their final location after the step. Use this when
  results storage is read-only on compute nodes. It is the default in
  `mfa_pipeline`.
* `execute_copy_inputs: True`: arubaito copies each step's inputs into `tmp_dir`
  before running it. Use this when the input data is not visible on the compute
  nodes at all.

With copying enabled, a placeholder like `{bam_file}` in a script refers to the
copy in `tmp_dir`, not to the final location. Scripts must therefore always use
the variables, never build file paths themselves. Later steps of the same job
reuse copies that already exist.

For local runs, set `execute_copy_outputs: False`, which saves the copying.
