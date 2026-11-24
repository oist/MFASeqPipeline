# Growth rate file format

An **input** file that the user provides, by default
`{meta_dir}/growthrates.tab` (variable `growthrates_file`). It gives the
exponential growth rate of each sample. The `velogp` and `mori` steps need it
(and the outdated `velocity_oscillation`).

Tab-separated text with no header. Each line has at least two columns:

| Column | Description |
|---|---|
| 1 | Sample name, exactly as the `sample` variable of the target |
| 2 | Growth rate λ in **1/hour** (population grows as `exp(λ t)`; doubling time = `ln 2 / λ`) |

Extra columns are ignored. Example:

```
exponential_sample_1	1.386
exponential_sample_2	0.693
```

The value is looked up when a step runs, through the variable
`growthrate_perhour`, which is defined as a shell command:
`$(awk -F'\t' '$1 == "{sample}" { print $2 }' '{growthrates_file}')`.
To skip the file, set `growthrate_perhour` directly, e.g. per target.
