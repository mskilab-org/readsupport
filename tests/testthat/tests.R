library(gUtils)
library(gGnome)
library(testthat)
library(Biostrings)
library(RSeqLib)


junction = jJ(parse.grl("1:1001-1001-,2:1001-1001-"))
reads = readRDS(system.file('extdata', 'reads.rds', package="readsupport"))
ref = readDNAStringSet(system.file('extdata', 'ref.fasta', package="readsupport"))
bwa = BWA(seq = ref)

test_that('junction.support', {

  res = junction.support(reads, junction, realign = FALSE)
  expect_equal(length(res), 49)
  expect_equal(length(unique(res$qname)), 18)

  ##  res = junction.support(reads, junction, ref = ref, bwa = bwa, realign = TRUE, pad.ref = 1e3)
  ##  expect_equal(length(res), 49)
  ##  expect_equal(length(unique(res$qname)), 18)
  
})

