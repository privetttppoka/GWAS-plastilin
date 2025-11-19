#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

df <- read.table(input_file, header=TRUE)

pvals <- df$P
pvals <- pvals[!is.na(pvals) & pvals > 0]

observed <- -log10(sort(pvals))
expected <- -log10(ppoints(length(pvals)))

png(output_file, width=1000, height=1000)

qqplot(expected, observed,
       xlab="Expected -log10(P)", ylab="Observed -log10(P)",
       main="QQ plot", pch=20)
abline(0, 1, col="red")

dev.off()
