#!/usr/bin/env Rscript

library(wasabi)

salmon_path <- '/stor/work/Lambowitz/cdw2854/bench_marking/alignment_free/salmon'
samples <- list.files(path = salmon_path, pattern = '[1-3]$')
sfdirs <- file.path(salmon_path, samples)
prepare_fish_for_sleuth(sfdirs)
