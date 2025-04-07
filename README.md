# MFASeqPipeline

`MFASeqPipeline` is a simple but flexiable pipeline for analyzing *marker frequency analysis* (MFA) data. The output of its first stage is a table of position-wise DNA abundances. From these abundances, the pipeline then computes the replication velocity as a function of location and time. The main features of `MFASeqPipeline` are

  - Supports non-parametric as well as parametric velocity estimation
  - Automatically detects and excludes unmappable regions
  - Works for arbitrary organisms as long as a reference genome is available

## Installation

## Quick-start guide

## The main pipeline steps

### _map_:  Mapping reads

### _fastqc_: Quality control

### _mapstats_: Mapping statistics

### _mappability_: Position-wise mappability

### _abundance_: Position-wise abundances

#### Output format

The output is a table containing the following columns:

1. *seq*: Chromosome (same name as in the reference genome for Yeast that we use)
2. *first*: Index of the first base in the window
3. *center*: Center of the window
4. *last*: Last base in the window
5. *windowsize*: Length of the window
6. *fragmentlength*: Fragment length of the sequenced fragments
7. *fragmentcoverage: Average number of bases of each fragment that overlap the current window (for small windows the window size, for large windows the fragment length, in between the value can be between these two values) sample_coverage_avg*: Average coverage of the bases in the window in the respective sample
8. *sample_coverage_sd*: Standard deviation of the coverage in the sample across the window
9. *sample_coverage_avg_stderr*: Estimated standard error of sample_coverage_avg
10. *sample_startingpoints*: Number of “starting points” in the sample, i.e. fragments whose first base is within the window
11. *sample_endpoints*: Number of “endpoints” in the sample, i.e. fragments whose last base is within the window
12. *sample_mappable*: 0 if there are regions in the sample within the window that are difficult to map, otherwise 1
13. *reference_coverage_avg*: Average coverage of the bases in the window in the respective reference example, i.e. in a non-growing population.
14. *reference_coverage_sd*: Standard deviation of the coverage in the reference example across the window
15. *reference_coverage_avg_stderr*: Estimated standard error of sample_coverage_avg.
16. *reference_startingpoints*: Number of “starting points” in the reference example, i.e. fragments whose first base is within the window
17. *reference_endpoints*: Number of “endpoints” in the reference example, i.e. fragments whose last base is within the window
18. *reference_mappable*: 0 if there are poorly mappable regions within the window, otherwise 1
19. *coverage_ratio*: Quotient sample_coverage_avg / reference_coverage_avg
20. *coverage_ratio_stderr*: Estimated standard error of coverage_ratio
21. *startingpoint_ratio*: Example starting points / reference starting points
22. *startingpoint_ratio_stderr*: Estimated standard error of startingpoint_ratio
23. *endpoint_ratio*: sample_endpoints / reference_endpoints
24. *endpoint_ratio_stderr*: Estimated standard error of endpoint_ratio.
25. *ratio_mappable*: 1 if sample_mappable == 1 and reference_mappable ==1, i.e. there are no mapping problems in either the sample or the reference sample within the window.

### __velocity_oscillation__: Frequency, phase and amplitude of velocity oscillations

### Handling single-end reads

## Additional modules

### _sra_: Downloading from the _Short Read Archive_ (SRA)

## The engine

### Variables

### Steps

### Execution

### Targets

