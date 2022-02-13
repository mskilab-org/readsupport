[![Build Status](https://travis-ci.org/mskilab/readsupport.svg?branch=master)](https://travis-ci.org/mskilab/readsupport)
[![Documentation Status](https://readthedocs.org/projects/gutils/badge/?version=latest)](https://readthedocs.org/projects/readsupport/?badge=latest)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/readsupport.svg)](https://codecov.io/github/mskilab/readsupport?branch=master)

readsupport
=======

Set of utility functions for use with `GenomicRanges`



Installation
------------

1. Install R-4.0 or later

2. Install devtools

```{r}
install.packages('devtools')
install.packages('testthat')
```
3. Install readsupport and dependent packages

```{r}
## allows dependencies that throw warnings to install
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)

devtools::install_github('mskilab/readsupport')
```

Read Support
-------------

What (short or long) reads support a given (known) rearrangement? This analysis is
important for vetting or merging rearrangement calls, for genotyping new
samples, and for co-calling.

`readsupport` allows you to assess the read support for a junction or arbitrary
rearranged contig (including those representing a multi-part rearrangements).
It is built on top of `gGnome`, `RSeqLib`, and `GenomicRanges`. 

Examples under construction!


Attributions
------------
> Marcin Imielinski - Associate Professor, Weill Cornell Medicine; Core Member,
> New York Genome Center

> Zi-Ning Choo - MD-PhD student, Tri-Institutional MD PhD program, Cornell
> University

> Alon Shaiber - Genomics Data Scientist, Imielinski Lab, Weill Cornell Medicine

[license]: https://github.com/mskilab/readsupport/blob/master/LICENSE


