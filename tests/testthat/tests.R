library(gUtils)
library(gGnome)
library(testthat)
library(Biostrings)
library(RSeqLib)

test_that('junction.support', {
    suppressWarnings(expr = {
        junction = jJ(parse.grl("1:1001-1001-,2:1001-1001-"))
        reads = readRDS(system.file('extdata', 'reads.rds', package="readsupport"))
        ref = readDNAStringSet(system.file('extdata', 'ref.fasta', package="readsupport"))
        bwa = BWA(seq = as.character(ref))

        res.norealn = junction.support(reads, junction, realign = FALSE, pad.ref = 1e3)
        expect_equal(length(res.norealn), 51)
        expect_equal(length(unique(res.norealn$qname)), 18)

        res = junction.support(reads, junction, ref = ref, bwa = bwa, realign = TRUE, pad.ref = 1e3)
        expect_true(length(res) > 30)
        expect_true(length(unique(res$qname)) > 15)
    })
})

## load additional unit testing cases
jj = gGnome::jJ(parse.grl("1:1001-1001+,2:1001-1001+"))

dup.fasta.fn = system.file("extdata", "dup", "ref.fasta", package = "readsupport")
dup.bwa = RSeqLib::BWA(seq = as.character(Biostrings::readDNAStringSet(dup.fasta.fn)))
dup.rs = Biostrings::readDNAStringSet(dup.fasta.fn)
dup.lifted.reads.gr = readRDS(system.file("extdata", "dup", "reads.rds", package = "readsupport"))

## positive control read qnames
dup.qn = c("C2UNVACXX131029:1:1313:13577:5314", "C2UNVACXX131029:1:2111:4328:55154", "C2UNVACXX131029:1:2315:7194:62620", "C2V4BACXX131029:8:1108:8496:70528", "C2V4BACXX131029:8:1110:9282:36756", "C2V4BACXX131029:8:1213:13865:72798", "C2V4BACXX131029:8:1313:3858:16480", "C2V4BACXX131029:8:2315:10689:35973", "C2V4CACXX131029:5:1104:19569:88882", "C2V4CACXX131029:5:1212:17513:49891", "C2V4CACXX131029:5:1304:10557:79471", "C2V4CACXX131029:5:1314:6774:11694", "C2V4CACXX131029:5:2114:12290:79915", "C2V4CACXX131029:5:2207:3605:84920", "C2V57ACXX131028:4:1115:13688:47062", "C2V57ACXX131028:4:1309:18079:19406", "C3416ACXX131026:3:2314:7168:82649")

test_that(desc = "test junction support for a lifted duplication", code = {
    suppressWarnings(
        expr = {
            dup.res = junction.support(reads = dup.lifted.reads.gr,
                                       junctions = jj,
                                       ref = dup.rs,
                                       bwa = dup.bwa,
                                       realign = TRUE,
                                       pad.ref = 1e3,
                                       verbose = TRUE)
            expect_true(length(intersect(dup.res$qname, dup.qn)) > 15)
        }
    )
})

tra.fasta.fn = system.file("extdata", "tra", "ref.fasta", package = "readsupport")
tra.bwa = RSeqLib::BWA(seq = as.character(Biostrings::readDNAStringSet(tra.fasta.fn)))
tra.rs = Biostrings::readDNAStringSet(tra.fasta.fn)
tra.lifted.reads.gr = readRDS(system.file("extdata", "tra", "reads.rds", package = "readsupport"))

tra.qn = c("C2UNVACXX131029:2:2212:5567:72624", "C2UNVACXX131029:2:2301:14123:12757", "C2UNVACXX131029:2:2314:8947:40374",  "C2UNVACXX131029:2:2315:16702:99138", "C2V4CACXX131029:6:1114:18447:13774", "C2V4CACXX131029:6:1114:7162:27490",  "C2V4CACXX131029:6:1207:10002:15275", "C2V4CACXX131029:6:1208:11658:27632", "C2V4CACXX131029:6:1307:18007:29051", "C2V4CACXX131029:6:1307:8739:73412",  "C2V4CACXX131029:6:1309:19784:74483", "C2V4CACXX131029:6:1309:4604:56089",  "C2V4CACXX131029:6:2112:18107:62803", "C2V4CACXX131029:6:2202:5390:9933",   "C2V4CACXX131029:6:2214:18785:78984", "C2V4CACXX131029:6:2310:4456:72900",  "C2V4CACXX131029:7:1106:14978:66273", "C2V4CACXX131029:7:1107:6694:56488",  "C2V4CACXX131029:7:1212:16437:70132", "C2V4CACXX131029:7:2106:11182:47280", "C2V4CACXX131029:7:2107:6308:43394", "C2V4CACXX131029:7:2109:18719:82584", "C2V4CACXX131029:7:2202:10962:74671", "C2V4CACXX131029:7:2202:5729:63480", "C2V4CACXX131029:7:2303:11673:17481", "C2V4CACXX131029:7:2311:18031:91318", "C2V57ACXX131028:3:1104:14587:32246", "C2V57ACXX131028:3:1109:16284:33989", "C2V57ACXX131028:3:1115:18086:95557", "C2V57ACXX131028:3:1301:13487:99357", "C2V57ACXX131028:3:2110:17865:44394", "C2V57ACXX131028:3:2111:14458:70988", "C2V57ACXX131028:3:2116:2239:19043",  "C2V57ACXX131028:3:2216:15163:17412", "C3416ACXX131026:4:1103:18472:60650", "C3416ACXX131026:4:1106:16255:68380", "C3416ACXX131026:4:1214:16987:30942", "C3416ACXX131026:4:2109:9907:38542",  "C3416ACXX131026:4:2116:19280:47641", "C3416ACXX131026:4:2305:8805:79929")

test_that(desc = "test junction support for a lifted translocation", code = {
    suppressWarnings(
        expr = {
            tra.res = junction.support(reads = tra.lifted.reads.gr,
                                       junctions = jj,
                                       ref = tra.rs,
                                       bwa = tra.bwa,
                                       realign = TRUE,
                                       pad.ref = 1e3,
                                       verbose = TRUE)
            expect_true(length(intersect(tra.res$qname, tra.qn)) > 15)
        }
    )
})

inv.fasta.fn = system.file("extdata", "inv", "ref.fasta", package = "readsupport")
inv.bwa = RSeqLib::BWA(seq = as.character(Biostrings::readDNAStringSet(inv.fasta.fn)))
inv.rs = Biostrings::readDNAStringSet(inv.fasta.fn)
inv.lifted.reads.gr = readRDS(system.file("extdata", "inv", "reads.rds", package = "readsupport"))

inv.qn = c( "HWI-D00222:200:C31LAACXX:2:1101:9642:48934", "HWI-D00222:200:C31LAACXX:2:1107:20393:80265", "HWI-D00222:200:C31LAACXX:2:1107:5592:96572",  "HWI-D00222:200:C31LAACXX:2:1115:12276:88474", "HWI-D00222:200:C31LAACXX:2:1203:13324:33883", "HWI-D00222:200:C31LAACXX:2:1214:15618:72475", "HWI-D00222:200:C31LAACXX:2:1304:17750:58748", "HWI-D00222:200:C31LAACXX:2:1310:2469:25826",  "HWI-D00222:200:C31LAACXX:2:1311:10933:17781", "HWI-D00222:200:C31LAACXX:2:1316:15402:56257", "HWI-D00222:200:C31LAACXX:2:2105:15113:23971", "HWI-D00222:200:C31LAACXX:2:2108:6826:97438",  "HWI-D00222:200:C31LAACXX:2:2109:2028:53934",  "HWI-D00222:200:C31LAACXX:2:2111:14886:78469", "HWI-D00222:200:C31LAACXX:2:2114:18158:3094",  "HWI-D00222:200:C31LAACXX:2:2115:5175:74274",  "HWI-D00222:200:C31LAACXX:2:2204:18359:27474", "HWI-D00222:200:C31LAACXX:2:2204:3031:96215",  "HWI-D00222:200:C31LAACXX:2:2207:6373:94363",  "HWI-D00222:200:C31LAACXX:2:2210:16836:77368", "HWI-D00222:200:C31LAACXX:2:2303:10186:91331", "HWI-D00222:200:C31LAACXX:2:2306:7093:33301",  "HWI-D00222:200:C31LAACXX:2:2315:9152:67614",  "HWI-ST1154:336:C31Y6ACXX:7:1101:12253:56264", "HWI-ST1154:336:C31Y6ACXX:7:1101:20935:48931", "HWI-ST1154:336:C31Y6ACXX:7:1108:19136:52786", "HWI-ST1154:336:C31Y6ACXX:7:1110:5823:33636",  "HWI-ST1154:336:C31Y6ACXX:7:1113:19197:70834", "HWI-ST1154:336:C31Y6ACXX:7:1202:11286:46381", "HWI-ST1154:336:C31Y6ACXX:7:1210:13004:72035", "HWI-ST1154:336:C31Y6ACXX:7:1214:19275:30134", "HWI-ST1154:336:C31Y6ACXX:7:1301:17873:30592", "HWI-ST1154:336:C31Y6ACXX:7:1310:3842:16324",  "HWI-ST1154:336:C31Y6ACXX:7:1310:5861:3838",   "HWI-ST1154:336:C31Y6ACXX:7:2104:5899:82406",  "HWI-ST1154:336:C31Y6ACXX:7:2115:10676:86109", "HWI-ST1154:336:C31Y6ACXX:7:2116:14221:14338", "HWI-ST1154:336:C31Y6ACXX:7:2116:14842:64468", "HWI-ST1154:336:C31Y6ACXX:7:2203:19847:43119", "HWI-ST1154:336:C31Y6ACXX:7:2210:20222:97211", "HWI-ST1154:336:C31Y6ACXX:7:2215:6363:74460",  "HWI-ST1154:336:C31Y6ACXX:7:2216:11146:8303",  "HWI-ST1154:336:C31Y6ACXX:7:2304:17698:92294", "HWI-ST1154:336:C31Y6ACXX:7:2312:5543:78833",  "HWI-ST1154:336:C31Y6ACXX:7:2316:20611:94788", "HWI-ST1353:296:C3391ACXX:5:1101:6039:52647",  "HWI-ST1353:296:C3391ACXX:5:1106:8511:10036",  "HWI-ST1353:296:C3391ACXX:5:1107:9173:53109",  "HWI-ST1353:296:C3391ACXX:5:1110:2036:52872",  "HWI-ST1353:296:C3391ACXX:5:1114:16055:86787", "HWI-ST1353:296:C3391ACXX:5:1201:17482:55204", "HWI-ST1353:296:C3391ACXX:5:1203:13645:50424", "HWI-ST1353:296:C3391ACXX:5:1205:13246:82380", "HWI-ST1353:296:C3391ACXX:5:1301:1445:46436",  "HWI-ST1353:296:C3391ACXX:5:1309:7584:45006",  "HWI-ST1353:296:C3391ACXX:5:1310:5200:11707",  "HWI-ST1353:296:C3391ACXX:5:2108:16137:26969", "HWI-ST1353:296:C3391ACXX:5:2113:3842:46874", "HWI-ST1353:296:C3391ACXX:5:2205:16594:47718", "HWI-ST1353:296:C3391ACXX:5:2209:16121:82636", "HWI-ST1353:296:C3391ACXX:5:2216:16415:64041", "HWI-ST1353:296:C3391ACXX:5:2306:17848:73247", "HWI-ST1353:296:C3391ACXX:5:2309:8861:79889",  "HWI-ST1353:296:C3391ACXX:5:2312:17808:44953", "HWI-ST898:499:C33AEACXX:2:1207:16865:11648",  "HWI-ST898:499:C33AEACXX:2:1216:3080:83398",   "HWI-ST898:499:C33AEACXX:2:1309:21304:51505",  "HWI-ST898:499:C33AEACXX:2:1316:13264:72133",  "HWI-ST898:499:C33AEACXX:2:2110:10998:13606",  "HWI-ST898:499:C33AEACXX:2:2201:13848:37777",  "HWI-ST898:499:C33AEACXX:2:2202:9232:5940",    "HWI-ST898:499:C33AEACXX:2:2203:3247:33503",   "HWI-ST898:499:C33AEACXX:2:2205:13495:76934",  "HWI-ST898:499:C33AEACXX:2:2207:18340:77518",  "HWI-ST898:499:C33AEACXX:2:2208:9012:9231",    "HWI-ST898:499:C33AEACXX:2:2210:7593:96885",   "HWI-ST898:499:C33AEACXX:2:2212:2641:25662",   "HWI-ST898:499:C33AEACXX:2:2216:17882:39040",  "HWI-ST898:499:C33AEACXX:2:2302:14479:4080",   "HWI-ST898:499:C33AEACXX:2:2302:8810:73137",   "HWI-ST898:499:C33AEACXX:2:2309:4182:10997", "HWI-ST898:499:C33AEACXX:2:2310:12494:94926", "HWI-ST898:499:C33AEACXX:2:2313:19134:48731")

test_that(desc = "test junction support for a lifted foldback-inversion", code = {
    suppressWarnings(
        expr = {
            inv.res = junction.support(reads = inv.lifted.reads.gr,
                                       junctions = jj,
                                       ref = inv.rs,
                                       bwa = inv.bwa,
                                       realign = TRUE,
                                       pad.ref = 1e3,
                                       verbose = TRUE)
            expect_true(length(intersect(inv.res$qname, inv.qn)) > 15)
        }
    )
})
