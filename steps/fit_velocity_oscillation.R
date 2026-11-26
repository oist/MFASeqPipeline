#!/usr/bin/env Rscript
library(data.table)
library(parallel)
library(this.path)
source(paste0(this.dir(), "/oscillating_velocity_models.R"))

GENOMES <- list(
  EColiMG1655=list(name="E. Coli MG1655", L=4641652, O=3925975)
)

MODELS <- list(
  vt=vt.harmonic.model,
  vd=vd.harmonic.model,
  vx=vx.harmonic.model
)

# Decode arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("fit_oscillation.R <genome> <model> <lambda> <abundance> <output> <cores>\n", file=stderr())
  quit(status=1)
}
argi <- 1
if (!(args[[argi]] %in% names(GENOMES)))
  stop("Unknown genome ", args[[argi]])
genome <- GENOMES[[args[[argi]]]]
argi <- argi + 1
if (!(args[[argi]] %in% names(MODELS)))
  stop("Unknown model ", args[[argi]])
model <- MODELS[[args[[argi]]]]
argi <- argi + 1
sample.lambda <- as.numeric(args[[argi]])
argi <- argi + 1
abundance.file <- args[[argi]]
if (!file.exists(abundance.file))
  stop("Abundance file ", abundance.file, " does not exist")
argi <- argi + 1
output.file <- args[[argi]]
if (file.exists(output.file))
  stop("Output file ", output.file, " already exist")
argi <- argi + 1
if (length(args) >= argi) {
  cores <- as.integer(args[[argi]])
} else {
  cores <- detectCores()
}
cat(paste0("Genome: ", genome$name, "\n"))
cat(paste0("Model: ", model$name, "\n"))
cat(paste0("Growth rate (1/h): ", sample.lambda, "\n"))
cat(paste0("Abundance: ", abundance.file, "\n"))
cat(paste0("Cores: ", cores, "\n"))

# Read abundances

d <- read.abundance(abundance.file, genome$L, genome$O)
cat(paste0("Datapoints: ", length(d$x), "\n"))

# Fit
ps <- fit.abundance.stepwise(d$x, d$y, d$w, L=genome$L, model=model,
                             lambda=sample.lambda/3600, mc.cores=cores)

# Build table of fitting results
tab <- rbindlist(lapply(names(ps), function(m) {
  p <- ps[[m]]
  df <- attr(p, "df")
  res <- attr(p, "residual")
  logl <- d$logl.const - res 
  as.data.table(c(list(model=model$abbreviation, regime=m,
                       n=length(d$x), df=df, optimal=FALSE,
                       L=genome$L, O=genome$O),
                  p[c("lambda", "D", "v0", "delta", "omega", "phi")],
                  list(residual=res, logl=logl, aic=2*df - 2*logl)))
}))

# Mark optimal model according to AIC
tab[order(aic, df)[1], optimal := TRUE]

# Output results
cat("Results:\n")
options(width=120)
print(tab, class=FALSE, row.names=FALSE, digits=3, scipen=0,
      nrows=nrow(tab))

# Write output
write.table(file=output.file, tab, sep="\t", quote=FALSE,
            col.names=TRUE, row.names=FALSE)
