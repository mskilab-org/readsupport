#' @import Biostrings
#' @import GenomicRanges
#' @import Matrix
#' @import data.table
#' @import gTrack
#' @import igraph
#' @import reshape2
#' @import bamUtils
#' @import gChain
#' @import gUtils
#' @import gGnome
#' @import RSeqLib
#' @import rtracklayer
#' @importFrom stats start end

#' @name junction.support
#' @title junction.support
#' @description
#'
#' Takes as input a GRanges of bam alignments (e.g. outputted from bamUtils::read.bam) and a GRanges of rearranged
#' reference aligned contigs (e.g. output of RSeqLib::BWA) and a set of Junction objects, and outputs reads supporting
#' these junctions by building a contig around each junction (from the reference) and then running contig.support (see
#' that functions docuemntation for criteria)
#'
#' @param reads GRanges in SAM / BAM format e.g. output of read.bam or BWA, with fields $qname, $cigar, $flag $seq all populated in standard fashion, and optionally $AS
#' @param junctions Junction object
#' @param bwa RSeqLib BWA object and path to fasta file corresponding to the reference
#' @param ref optional DNAStringSet corresponding to reference genome sequence
#' @param pad padding around the junction breakpoint around  which to analyze contig and reference sequences, this should be several standard deviations above the average insert size (2000)
#' @param realign flag whether to realign or just use existing alignments
#' @param bx logical flag whether data is linked reads, must then have BX flag, and the pad will be set to minimum 1e5
#' @param pad.ref (numeric) pad around breakend for grabbing reference sequence, default pad * 20
#' @param both (logical) perform both realignment and non-realignment-based support? (default TRUE)
#' @param walks (gWalk) walks corresponding to junction (default NULL)
#' @param verbose logical flag (TRUE)
#' @param ... additional parameters to contig support
#' @return reads re-aligned to the reference through the contigs with additional metadata describing features of the alignment
#' @export
#' @author Marcin Imielinski
junction.support = function(reads,
                            junctions = NULL,
                            bwa = NULL,
                            ref = NULL,
                            pad = 500,
                            bx = FALSE,
                            pad.ref = pad*20,
                            both = TRUE,
                            realign = TRUE,
                            walks = NULL,
                            verbose = TRUE, ...)
{

  if (!inherits(reads, 'GRanges') || is.null(reads$qname) || is.null(reads$cigar) || is.null(reads$seq) || is.null(reads$flag))
    stop('read input must be GRanges with fields $qname, $cigar, $seq, $flag and optionally $AS')

  if (bx)
    pad = max(pad, 1e5)

  if (!is.null(junctions))
    walks = gGnome::jJ(junctions$grl)$gw(pad = pad)

  if (is.null(walks))
    stop('Either walks or junctions must be provided')

  if (bx)
  {
    if (is.null(reads$BX))
      stop('reads must have BX tag, may need to read.bam with tag option to extract it')

    if (!length(reads))
      return(reads)

    sc = score.walks(walks$grl, reads = reads, verbose = FALSE, raw = TRUE)$sc
    res = as.data.table(melt(as.matrix(sc)))[value>0, .(BX = Var1, walk = Var2)]
    reads = gr2dt(reads) %>% merge(res, by = 'BX') %>% dt2gr
    return(reads)
  }

  if (!realign)
  {
    if (is.null(junctions))
      junctions = walks$edges$junctions

    ## strand flip since 
    ## read orientation convention
    ## is opposite to junction convention
    reads = gr.flipstrand(reads) 
    reads$R1 = bamUtils::bamflag(reads$flag)[,'isFirstMateRead']>0
    r1 = reads %Q% (R1 == TRUE) %>% as.data.table
    r2 = reads %Q% (R1 == FALSE) %>% as.data.table
    ov = merge.data.table(r1, r2, by = 'qname')
    if (!nrow(ov))
      return(reads[c()])
    
    sl = seqlengths(reads)
    grl = grl.pivot(
      GRangesList(dt2gr(ov[, .(seqnames = seqnames.x, start = start.x, end =end.x, strand = strand.x)],
                        seqlengths = sl),
                  dt2gr(ov[, .(seqnames = seqnames.y, start = start.y, end = end.y, strand = strand.y)],
                        seqlengths = sl)))
    values(grl)$qname = ov$qname
    ## make junctions out of reads and cross with "real" junctions
    jn = gGnome::merge(jJ(grl), junctions, cartesian = TRUE, pad = pad)
    if (!length(jn))
      return(reads[c()])
    out = merge(as.data.table(gr.flipstrand(reads)), unique(jn$dt[, .(qname, junction.id = subject.id)]), by = 'qname') %>% dt2gr(seqlengths = sl)
    return(out)
  }

  if (is.null(bwa) & is.null(ref))
    stop('BWA object or reference must be provided if realign = TRUE')

  if (is.null(bwa) && !is.null(ref) && is.character(ref))
    bwa = BWA(ref)

  if (is.null(bwa))
    stop('BWA object or reference must be provided if realign = TRUE')
  
  if (inherits(bwa, 'character') && file.exists(bwa))
  {
    if (verbose)
      message('Loading BWA index')
    bwa = BWA(bwa)
  }

  if (!inherits(ref, 'DNAStringSet'))
  {
    if (verbose)
      message('Loading genome reference as DNAStringSet')

    ref = rtracklayer::import(bwa@reference)
  }

  ## only >use the fasta header before the first space as the seqnames of ref 
  names(ref) = strsplit(names(ref), '\\s+') %>% sapply('[', 1)

  if (length(setdiff(seqnames(walks$nodes$gr), seqlevels(ref))))
    stop('seqlevels mismatch between junctions / walks and reference, please check reference (e.g. chr issues)')

  if (length(setdiff(seqnames(walks$nodes$gr), seqlevels(bwa))))
    stop('seqlevels mismatch between junctions / walks and BWA reference, please check reference (e.g. chr issues)')

  if (verbose)
    message('Building and mapping derivative contigs')

  contig = bwa[ref[gr.fix(walks$grl, ref, drop = TRUE)]]

  if (verbose)
    message('Building reference contigs flanking junctions')
  contigref = ref[gr.fix(walks$edges$junctions$footprint + pad.ref, ref, drop = TRUE)]


  if (verbose)
    message('Making gChain mapping contigs to reference')
  cg.contig = gChain::cgChain(contig)

  if (verbose)
    message('Running contig support')

  reads = contig.support(reads, contig, ref = bwa, cg.contig = cg.contig, ...)
  ##  reads = contig.support(reads, contig, ref = contigref, cg.contig = cg.contig, ...)
  reads$junction.id = as.integer(as.character(reads$contig.id))

  if (both)
  {
    if (length(reads))
      reads$source = 'realignment'
    
    reads2 = junction.support(reads, junctions, bwa = bwa, ref = ref, pad = pad, realign = FALSE, both = FALSE, pad.ref = pad.ref, walks = walks, verbose = verbose, ...)

    if (length(reads2))
      reads2$source = 'original_alignment'

    reads = grbind(reads, reads2)
  }
  return(reads)  
}


#' @name contig.support
#' @title contig.support
#' @description
#'
#' Takes as input a GRanges of bam alignments (e.g. outputted from bamUtils::read.bam) and a GRanges of rearranged
#' reference aligned contigs (e.g. output of RSeqLib::BWA).
#'
#' It identifies the subset of reads that support each of the contigs and "lifts" those reads
#' through the read --> contig and contig --> reference alignments, returning supporting reads in reference coordinates.
#'
#' The criteria for support include min.bases aligning to at least two chunks of the rearranged contig, and
#' requirement that min.aligned.frac fraction of bases in every supporting read is aligned to that contig.
#'
#' Additional requirements for support include not allowing split alignment of individual reads to the contigs
#' (note: this does not mean we don't detect split reads that support the structural variant, this is captured
#' by the contig -> reference alignment, we are just requiring the reads align (near) perfectly to the contig).
#' and requiring alla alignments from a read pair (oriented to R1 frame of the fragment) to align to the same
#' strand of the contig.
#'
#' Finally, reads are not included in support if they align better to the reference than their native alignment,
#' which is determined by comparing the $AS of their contig alignment with their original alignment score, stored
#' in the provided metadata $AS field.  If reference AS is not provided as metadata, it will is assumed to be zero. 
#'
#' $AS can be optionally recomputed against a DNAStringSet "ref" that represent the reference
#' sequence.  (Note that this "ref" does not have to be the full genome reference, it is just used to compute
#' the alignment scores, and in fact for this to work  efficiently, it's recommended that the provided
#' reference sequence is local to the regions of interest, e.g. a few kb flanking each SV breakpoint,
#' rather than the whole genome.)
#'
#' The outputted reads include additional metadata including number of bases aligning to each chunk of the aligned contig.
#' 
#' 
#' @param reads GRanges in SAM / BAM format e.g. output of read.bam or BWA, with fields $qname, $cigar, $flag $seq all populated in standard fashion, and optionally $AS
#' @param contig GRanges in SAM / BAM format wth fields $qname, $cigar and $seq all [populated
#' @param ref optional DNAStringSet representing a reference sequence to compute alignments against
#' @param chimeric logical flag whether to require reads to support junctions in chimericcontigs (ie discontiguous chunks on the reference), chimeric = FALSE
#' @param strict strict requires that the alignment score of the read to contig alignment needs to be better for at least one read (and also not worse for any of the reads)
#' @param cg.contig (cgChain) gChain corresponing to contig (maps contig coordinates to ref coordinates)
#' @param isize.diff (numeric) the insert size in the reference vs. the insert size in the contig between discordant read pairs
#' @param min.bases (numeric) min aligned bases to be considered valid alignment, default 20
#' @param min.aligned.frac (numeric) min fraction of bases in alignment that are matching, default 0.95
#' @param new (logical) new scoring scheme, default TRUE
#' @param bowtie (logical) use bowtie for read alignment? requires bowtie to be callable from command line, default FALSE (use BWA)
#' @param outdir (logical) output directory for bowtie2 temporary files
#' @param verbose (logical) print stuff (default TRUE)
#' 
#' @return reads re-aligned to the reference through the contigs with additional metadata describing features of the alignment
#' @export
#' @author Marcin Imielinski

contig.support = function(reads,
                          contig,
                          ref = NULL,
                          chimeric = TRUE,
                          strict = TRUE,
                          cg.contig = gChain::cgChain(contig),
                          isize.diff = 1e3,
                          min.bases = 20,
                          min.aligned.frac = 0.95,
                          new = TRUE,
                          bowtie = FALSE,
                          outdir = "./",
                          verbose = TRUE)
{
    if (length(reads)==0)
        stop('reads must be non empty GRanges with $qname, $cigar, $seq, and $flag fields')

    if (length(contig)==0)
        stop('contigs must be non empty GRanges with $qname, $cigar and $seq fields')

    if (verbose)
        message('Prepping reads for contig alignment')
    seq = unique(gr2dt(contig), by = c('qname'))[, structure(as.character(seq), names = as.character(qname))]
    bwa.contig = RSeqLib::BWA(seq = seq)
    chunks = gChain::links(cg.contig)$x
    strand(chunks) = '+'
    chunks = disjoin(chunks)
    if (is.null(reads$R1))
        reads$R1 = bamUtils::bamflag(reads$flag)[,'isFirstMateRead']>0
    reads$read.id = 1:length(reads)
    if (is.null(reads$AS))
    {
        warning('AS not provided in reads, may want to consider using tag = "AS" argument to read.bam or provide a ref sequence to provide additional specificity to the contig support')
        reads$AS = 0
    }
    reads = reads %Q% (!is.na(seq))##[!is.na(seq)]
    nix = as.logical(strand(reads) == '-' )
    reads$seq[nix] = reverseComplement(DNAStringSet(reads$seq[nix])) ## flip read sequences to original strand
    reads[!reads$R1] = gr.flipstrand(reads[!reads$R1]) ## flip R2 read orientation to R1 strand
    reads$seq[!reads$R1] = reverseComplement(DNAStringSet(reads$seq[!reads$R1])) ## flip R2 read sequences to R1 strand
    ## remove non-primary alignments?
    ##reads = reads %Q% (bamUtils::bamflag(reads$flag)[,'isNotPrimaryRead']==0)
    reads = reads %Q% (!duplicated(paste(qname, R1))) ## this will remove split reads??

    if (!is.null(ref)) ## realign reads against reference DNAStringSet if provided to get alignment scores
    {
        if (verbose)
            message('Realigning reads against reference DNAStringSet')
        
        if (inherits(ref, 'character') | inherits(ref, 'DNAStringSet')) {
            if (verbose) {
                message("building BWA from sequence")
            }
            bwa.ref = RSeqLib::BWA(seq = as.character(ref))
        } else if (inherits(ref, 'BWA')) {
            if (verbose) {
                message("Using supplied BWA")
            }
            bwa.ref = ref
        } else {
            stop("invalid data type provided for ref: ", class(ref)[1])
        }

        ## tmp = bwa.ref[reads$seq] %>% gr2dt
        if (bowtie) {
            if (verbose) { message("Using Bowtie for read to contig alignment!") }
            ## infer the name from the basename of the ref file
            if (!inherits(ref, 'character')) {
                stop("Must supply path to reference")
            }
            fasta.dirname = dirname(ref)
            fasta.basename = basename(ref)
            if (!dir.exists(fasta.dirname)) {
                stop("Ref directory does not exist!")
            }
            fasta.name = gsub("\\.(fasta|fa|fastq|fq).*$", "", fasta.basename)
            if (verbose) { message("Using reference name: ", fasta.name) }
            tmp = bowtie_aln(reads = as.data.table(reads),
                             ref.dir = dirname(ref),
                             ref.basename = fasta.name,
                             outdir = normalizePath(file.path(outdir, "ref")),
                             verbose = verbose)
            tmp = as.data.table(tmp)
        } else {
            tmp = bwa.ref[reads$seq] %>% gr2dt
            tmp$ix = as.numeric(as.character(tmp$qname))
            tmp$R1 = reads$R1[tmp$ix]
            tmp$qname = reads$qname[tmp$ix]
        }
        tmp = unique(tmp, by = c('qname', 'R1'))
        setkeyv(tmp, c('qname', 'R1'))
        if (nrow(tmp))
        {
            ## what if not both read pairs are included?
            tmp[, isize := ifelse(any(seqnames != seqnames[1] | any(strand != strand[1])), Inf, diff(range(start, end))), by = qname]
            ## isize should be NA if only one of the sides has a match
            unpaired.qnames = tmp[, .(count = .N), by = qname][count == 1, qname]
            tmp[qname %in% unpaired.qnames, isize := NA_integer_]
            ## the final isize should only account for the isize from the provided refseq
            reads$isize = tmp[.(reads$qname, reads$R1), isize]## pmin(tmp[.(reads$qname, reads$R1), isize], Inf, na.rm = TRUE)
            ## next line is the original version
            ## reads$isize = pmin(tmp[.(reads$qname, reads$R1), isize], reads$isize, Inf, na.rm = TRUE)
            ## reads[, c("isize", "isize.new", "isize.pmin", "sample")] %Q% (sample == "A0K8N")
            ## reads %Q% (qname %in% unpaired.qnames)
            reads$AS = tmp[.(reads$qname, reads$R1), AS]
            ## replace cigar
            reads$ref.aligned.cigar = tmp[.(reads$qname, reads$R1), cigar]
        }
    }

    if (verbose)
        message('Aligning reads against derivative contigs')
    

    ## aligning reads to contig
    rdt = as.data.table(reads)
    ## use new cigar?
    rdt[, ref.aligned := countCigar(ref.aligned.cigar)[, 'M']] 
    rdt[, ref.aligned := countCigar(cigar)[, 'M']] ## this is not the cigar from the realignment?
    rdt[, ref.aligned.frac := ref.aligned/qwidth[1], by = .(qname, R1)]

    reads$ref.aligned.frac = rdt$ref.aligned.frac

    ## readsc is a data.table that denotes read locations in contig coordinates
    if (bowtie) {
        contigs.dt = unique(gr2dt(contig), by = c('qname'))
        readsc = bowtie_aln(reads = as.data.table(reads),
                            ref.seq = contigs.dt[, .(qname, seq)],
                            outdir = normalizePath(file.path(outdir, "contig")),
                            verbose = verbose)
        readsc = as.data.table(readsc)
    } else {
        readsc = bwa.contig[reads$seq] %>% gr2dt
        readsc$ix = as.integer(as.character(readsc$qname))
        readsc$R1 = reads$R1[readsc$ix]
        
    }

    ## important! if no reads align to the contig then the following will fail!
    if (!nrow(readsc)) {
        return(reads[c()])
    }

    readsc$read.id = reads$read.id[readsc$ix]
    readsc$cigar = as.character(readsc$cigar)
    

    ## these are splits on the contig, not reference --> shouldn't be any for good alignment
    ## readsc[, nsplit := .N, by = .(qname, R1)] 
    readsc[, aligned := countCigar(cigar)[, 'M']]

    ## these are splits on the contig, not reference --> shouldn't be any for good alignment
    readsc[, aligned.frac := aligned/qwidth[1], by = .(qname, R1)]
    readsc$AS.og = reads$AS[readsc$ix]
    readsc$isize = abs(reads$isize[readsc$ix])


    readsc$seqnames.og = seqnames(reads)[readsc$ix] %>% as.character
    readsc$strand.og = strand(reads)[readsc$ix] %>% as.character
    readsc$start.og = start(reads)[readsc$ix]
    readsc$end.og = end(reads)[readsc$ix]
    readsc$ref.isize = reads$isize[readsc$ix]
    ## readsc$ref.isize = gr2dt(readsc)[, ref.isize := ifelse(
    ##                                        both &
    ##                                        all(seqnames.og == seqnames.og[1]) &
    ##                                        all(strand.og == strand.og[1]),
    ##                                        as.numeric(diff(range(c(start.og, end.og)))),                                   
    ##                                        Inf), by = qname]$ref.isize %>% abs
    readsc$ref.aligned.frac = reads$ref.aligned.frac[readsc$ix]
    readsc$AS.og[is.na(readsc$AS.og)] = 0
    readsc$qname = reads$qname[readsc$ix]

    ## track sample (comment out later)
    ## tst = readsc[, .(aligned.frac, AS, AS.og, ref.aligned.frac, qname)][, sample := rdt$sample[match(qname, rdt$qname)]]
    ## tst[, .(aligned.frac, ref.aligned.frac, AS, AS.og, sample)]
    
    ## new scoring method based on cgChain of reads to contigs
    if (new)
    {
        ## cgChain representing read to contig alignments
        readsc$al.id = 1:nrow(readsc)

        if (verbose)
            message('Generating read to contig cgChain')
        alcg = gChain::cgChain(readsc)
        alchunks = cbind(as.data.table(values(alcg)),
                         as.data.table(gChain::links(alcg)$x),
                         as.data.table(gChain::links(alcg)$y)[, .(contig = seqnames,
                                                                  contig.start = start,
                                                                  contig.end = end,
                                                                  contig.strand = strand)])

        ## strands should be aligned to read / fragment + strand, but if not let's flip
        alchunks[strand == '-', ":="(strand = '+', contig.strand = c('+' = '-', '-' = '+')[contig.strand])]

        ## now for each al.id (ie bam record) let's pick the left most gChain / links record on the read / fragment 
        ## ie this is the lowest coordinate on the query
        ## (note that cgChain will split indels into separate ranges hence giving one to many mapping of al.id
        ## to records in links)
        ## setkeyv(alchunks, c('qname', 'start', 'end'))
        setkeyv(alchunks, c('contig', 'contig.start', 'contig.end'))
        ## alchunks[, is.min := start == min(start), by = al.id]
        ## alchunks = alchunks[is.min == TRUE, ]

        ## so now we want to find alignments that are
        ## (1) concordant with respect to the contig
        ##  i.e. there is a monotonic increase (decrease) of contig.start if the contig.strand is + (-)
        ## (2) most of the read (aligned.frac) is represented
        ## (3) AS scores are better than original
        ## (4) isize better than original (where isize is the contig. span between the first and last alignment) .. related to (1)


        if (verbose)
            message('Scoring read to contig to alignments')
        alchunks[, contig.sign := ifelse(contig.strand == '+', 1, -1)]
        alchunks[, concordant.sign := all(contig.sign == contig.sign[1]), by = qname]

        ## check if we never go from R1 == FALSE to R1 == TRUE
        ## alchunks[, concordant.R1R2 := all(diff(!R1)>=0), by = qname]

        ## check to see that our contig.start always increasing or decreasing
        ## alchunks[, concordant.start := all(diff(contig.sign[1]*contig.start)>0), by = qname]

        ## new code from marcin to make concordant.R1R2 and concordant.start robust to small inserts
        alchunks[, both := any(R1) & any(!R1), by = qname]
        alchunks[, concordant.R1R2 := ifelse(both,(contig.sign*sign((contig.start[R1][1]<contig.start[!R1][1]) - 0.5))>0, TRUE), by = qname]
        ##alchunks[, concordant.start := all((contig.sign[1]*diff(start))>0), by = .(qname, R1)]
        alchunks[, concordant.start := all((contig.sign[1]*diff(start))>0), by = .(qname, R1)]

        ## alchunks[, contig.isize := ifelse(both,diff(range(contig.start, contig.end)), 0), by = qname]
        alchunks[, contig.isize := ifelse(all(contig == contig[1]) &
                                          all(contig.strand == contig.strand[1]),
                                          abs(as.numeric(diff(range(c(contig.start, contig.end))))),
                                          Inf), by = qname]
        alchunks[(!both), contig.isize := NA_integer_]
        alchunks[, bases := sum(width), by = .(qname)]

        ## nonzero width of chunks with better alignment scores relative to REF
        ## should this also be by read? (e.g. qname, R1)
        alchunks[, AS.better := sum(width[AS>AS.og]), by = .(qname)]
        ## keep qname if a single read is truncated
        alchunks[, AS.worse := sum(width[AS<AS.og]), by = .(qname)] 
        alchunks[, AS.equal := sum(width[AS==AS.og]), by = .(qname)]
        alchunks[, nsplit := .N, by = .(qname, R1)]

        keepq = alchunks[concordant.sign & concordant.R1R2 & concordant.start &
                         bases > min.bases & aligned.frac > min.aligned.frac &
                         aligned.frac >= ref.aligned.frac &
                         nsplit  == 1 & ## no split reads
                         ((contig.isize - ref.isize < isize.diff) |
                          (is.infinite(ref.isize) & !is.infinite(contig.isize)) |
                          (is.na(ref.isize) & !is.infinite(contig.isize))) & 
                         (AS.better>0 | (((ref.isize > 1e3 & ref.isize > contig.isize) | is.infinite(ref.isize) | is.na(ref.isize)) & (!is.infinite(contig.isize) & !is.na(contig.isize)))) & 
                         AS.worse == 0, ]

        ## keep read-specific information...
        keepq = keepq[, .(qname, contig, contig.id = as.character(contig), contig.isize, contig.strand, bases, contig.sign, AS.better, AS.worse, AS.equal)] %>% unique(by = 'qname')
        ## keepq = keepq[, .(qname, R1, contig, contig.isize, contig.strand, bases,
        ##                   contig.sign, AS.better, AS.worse, AS.equal)] ## %>% unique(by = c('qname', 'R1'))
    }
    else ## old scoring method
    {
        ## if strict (default) remove any alignments that overlap others in the same qname
        if (strict)
        {
            readsc = dt2gr(readsc)
            readsc = readsc %Q% (rev(order(AS)))
            readsc = readsc[!duplicated(gr.match(readsc, readsc, by = 'qname')), ] %>% gr2dt
        }


        if (verbose)
            message('Computing overlap stats')

        ov = dt2gr(readsc) %*% chunks
        strand(ov) = readsc$strand[ov$query.id]
        ov$subject.id = paste0('chunk', ov$subject.id)
        ovagg = dcast.data.table(ov %>% gr2dt, qname ~ subject.id, value.var = 'width', fun.aggregate = sum)
        ovagg$nchunks = rowSums(ovagg[, -1]>min.bases)  ## good means we hit multiple chunks with sufficient bases
        rstats = gr2dt(ov)[, .(
                          contig.id = unique(seqnames)[1],
                          pos = sum(width[strand == '+']),
                          neg = sum(width[strand == '-']),
                          aligned.frac = min(aligned.frac),
                          num.contigs = length(unique(seqnames)), ### fixing later ... multiple contigs as input could distort results
                          paired = any(R1) & any(!R1), 
                          isize.contig = diff(range(c(start, end))),
                          isize.og = isize[1],
                          qsplit = any(nsplit>1), ## any sequences in this qname split on the contig ie a bad alignment on the contig
                          worse = any(AS.og>AS), ## any alignment in this qname worse than vs reference?
                          better = any(AS>AS.og) ## any alignment in this qname better than vs reference?
                      ), by = qname] %>% merge(ovagg, by = 'qname')

        ## apply filters ie nchunks>1 if chimeric, all alignments have to be of one sign
        ## if not paired then AS < AS.og else isize<isize.og
        keepq = rstats[nchunks>chimeric & (pos == 0 | neg  == 0) & aligned.frac > min.aligned.frac & !worse & (better | !strict | (paired & isize.contig < isize.og - isize.diff)) & !qsplit & num.contigs == 1, ]
        if (nrow(keepq)==0)
            return(reads[c()])

        keepq$aligned.frac = NULL
    }

    ## R1/R2-speicific information is lost here
    readsc = merge(readsc, keepq, by = 'qname') %>% dt2gr
    ## readsc = merge(readsc, keepq, by = c('qname', 'R1')) %>% dt2gr
    
    if (verbose)
        message('Lifting reads through contig back to reference')

    out = gChain::lift(cg.contig, readsc)

    if (length(out)) ## add reads metadata back to out
    {
        out[!out$R1] = gr.flipstrand(out[!out$R1])
        out$col = ifelse(out$R1, 'blue', 'gray')

        if (verbose)
            message('Adding metadata to reads')
        metacols = setdiff(names(values(reads)), names(values(out)))
        values(out) = cbind(values(out), values(reads)[match(out$read.id, reads$read.id), metacols])
    }

    if (verbose)
        message('Done')

    return(out)
}
