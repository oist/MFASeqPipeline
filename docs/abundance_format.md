# Abundance file format

Produced by the `abundance` step of [mfa_pipeline](mfa_pipeline.md) and
[mfa_se_pipeline](mfa_se_pipeline.md) (`steps/abundance.py`). It is the input
for `velogp`, `abundance_plot` and the outdated `velocity_oscillation`.

* Tab-separated text, gzip-compressed if the file name ends in `.gz` (the
  default name is `{out_dir}/{sample}.abd.w{windowsize}bp.gz`).
* The first line is a header with the column names below.
* One row per window. Each chromosome is split into consecutive,
  non-overlapping windows of `windowsize` bases. The last window of a
  chromosome may be shorter.
* Coordinates are **0-based**, and `first`/`last` are **inclusive**.
* Rows are grouped by chromosome, in the order of the BAM header.

Only fragments with mapping quality ≥ `mapping_quality` are counted. In
paired-end mode a fragment covers the bases from the start of the forward read
to the end of the reverse read. In single-end mode only the read itself counts.

## Columns

| # | Column | Type | Description |
|---|---|---|---|
| 1 | `seq` | string | Chromosome name, as in the reference genome |
| 2 | `first` | int | First base of the window (0-based) |
| 3 | `midpoint` | float | Center of the window, `(first + last) / 2` |
| 4 | `last` | int | Last base of the window (0-based, inclusive) |
| 5 | `windowsize` | int | Length of the window in bases |
| 6 | `fragmentlength` | float | Mean fragment length in the sample, for this chromosome |
| 7 | `fragmentcoverage` | float | Expected number of bases of a fragment that overlap the window, given that it overlaps at all. This is about the window size for small windows and the fragment length for large windows |
| 8 | `sample_coverage_avg` | float | Mean per-base coverage of the window in the sample |
| 9 | `sample_coverage_sd` | float | Standard deviation of the per-base coverage within the window |
| 10 | `sample_coverage_avg_stderr` | float | Estimated standard error of `sample_coverage_avg`, assuming Poisson-distributed fragment counts: `sqrt(avg) * sqrt(fragmentcoverage / windowsize)` |
| 11 | `sample_startingpoints` | int | Number of fragments whose first base lies in the window |
| 12 | `sample_endpoints` | int | Number of fragments whose last base lies in the window |
| 13 | `sample_mappable` | 0/1 | 1 if the window does not overlap any non-mappable region of the sample's [mappability file](mappability_format.md) |
| 14–19 | `reference_*` | | Same as columns 8–13, for the reference (non-replicating) sample |
| 20 | `coverage_ratio` | float | `sample_coverage_avg / reference_coverage_avg`: the relative DNA abundance |
| 21 | `coverage_ratio_stderr` | float | Estimated standard error of `coverage_ratio` (see below) |
| 22 | `startingpoint_ratio` | float | `sample_startingpoints / reference_startingpoints` |
| 23 | `startingpoint_ratio_stderr` | float | Estimated standard error of `startingpoint_ratio` |
| 24 | `endpoint_ratio` | float | `sample_endpoints / reference_endpoints` |
| 25 | `endpoint_ratio_stderr` | float | Estimated standard error of `endpoint_ratio` |
| 26 | `ratio_mappable` | 0/1 | `sample_mappable AND reference_mappable` |

The standard errors of the ratios use the following approximation for the
coefficient of variation `CV = stderr / mean` of a quotient `X / Y`:

```
CV[X/Y]^2 = (CV[X]^2 + CV[Y]^2 - CV[X]^2 * CV[Y]^2) / (1 - 2 * CV[Y]^2)
```

For start and end points, Poisson counts are assumed, so `CV^2 = 1 / count`.

## Without a reference sample

If no reference is given (e.g. with `mfa_no_reference.yaml`), the reference and
ratio columns are still present but hold placeholders. The `reference_coverage_*`
columns and all ratio columns are `nan`, `reference_startingpoints` and
`reference_endpoints` are `-1`, and `reference_mappable` and `ratio_mappable` are
`0`. Downstream steps must then use the `sample_*` columns, which
`mfa_no_reference.yaml` arranges for `velogp` and `abundance_plot`.

## Example

```
seq	first	midpoint	last	windowsize	fragmentlength	fragmentcoverage	sample_coverage_avg	...
NC_000913.3	0	4999.5	9999	10000	312.00	311.99	153.21	...
NC_000913.3	10000	14999.5	19999	10000	312.00	311.99	149.87	...
```
