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

  res.norealn = junction.support(reads, junction, realign = FALSE, pad.ref = 1e3)
  expect_equal(length(res.norealn), 51)
  expect_equal(length(unique(res.norealn$qname)), 18)

  res = junction.support(reads, junction, ref = ref, bwa = bwa, realign = TRUE, pad.ref = 1e3)
  expect_true(length(res) > 30)
  expect_true(length(unique(res$qname)) > 15)
})

## sdqn = setdiff(res.norealn$qname, res$qname)
## sdres = res.norealn %Q% (qname %in% sdqn)

## grl = split(res, res$qname)
## sdgrl = split(sdres, sdres$qname)

## rg = stack(junction$grl)
## gw = gW(grl = grl)
## sdgw = gW(grl = sdgrl)
## gt = gTrack(grl, labels.suppress.grl = TRUE)
## sdgt = gTrack(sdgrl, labels.suppress.grl = TRUE)

## ppng(plot(c(sdgt, gt), rg + 1e3, links = junction$grl))


  ##  expect_equal(length(res), 49)
  ##  expect_equal(length(unique(res$qname)), 18)
  
## })

