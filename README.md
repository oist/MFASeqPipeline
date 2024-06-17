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

### __velocity_oscillation__: Frequency, phase and amplitude of velocity oscillations

### Handling single-end reads

## Additional modules

### _sra_: Downloading from the _Short Read Archive_ (SRA)

## The engine

### Variables

### Steps

### Execution

### Targets

