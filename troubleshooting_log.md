# Log 10/04/23

The major problem that exists with the Parsebio run is that barcode one is matching very few reads. It is the first read to be added and therefore unlikely to be added at low efficiency. Furthermore it has a linker right before it which can be used to calibrate the read position even if some bases were skipped.\
\
In order to better understand the makeup of `read 2`, which contains the barcode and UMI reads, we are going to load the data into R. Will try to answer the following questions with our analyses.

- [x] What is the distribution of read lengths?
- [x] What percent of reads have `linker B` and and `linker A` in the correct position?
- [x] What are the barcodes when using the linkers as strict guides to anchor barcode extraction?
    - cross reference this to list of barcodes to see is the right ones are really being represented.
- [x] What are the barcodes when using the presupposed character indices to determine where the barcodes should be?
    - cross reference this to list of barcodes to see is the right ones are really being represented.

Going to move forward with `pull-sublibrary-5_S3_R2_001.fastq.gz`.

```{.r}
# use microseq package to read in file
renv::install('microseq') # for fastq
renv::install('R.utils') # for gzipped files
# read in file
library('microseq')
WD = '/home/patrick/test/01-basecalled/'
file = 'pull-sublibrary-5_S3_R2_001.fastq.gz'
fq = readFastq(paste0(WD, file))
# all sequences are contained in fq$Seqeuence. 178162 sequences
# determine length distribution of reads
fqL = table(unlist(lapply(fq$Sequence, nchar)))
fqL = 100*fqL/sum(fqL)
fqL = fqL[-which(fqL<1)]
#       83        84        85        86 
# 2.237851  3.627597 27.630471 64.174740
# 92% of reads are 85 or 86 bp long
# What percent of reads have `linker B` and and `linker A` in the correct position?
linkerB = lapply(fq$Sequence, substr, 19, 48)
100*table(unlist(lapply(linkerB, function(x) sum(strsplit(x, '')[[1]] != strsplit('GTGGCCGATGTTTCGCATCGGCGTACGACT', '')[[1]]))))/length(fq$Sequence)
#           0           1           2           3           4           5           6           7           8           9 
# 80.48012483  2.60830031  0.33284314  0.24528238  0.22675992  0.14312816  0.18410211  0.18971498  0.18690854  0.26661129 
#          10          11          12          13          14          15          16          17          18          19 
#  0.46081656  0.63088650  0.71395696  0.52704842  0.42433291  0.30870781  0.34126245  0.44117152  0.59384156  0.69262806 
#          20          21          22          23          24          25          26          27          28          29 
#  1.00751002  1.31060496  1.52670042  1.61931276  1.45429441  1.82811149  0.60955759  0.28962405  0.17792795  0.10215422 
# 80% have linker B in the right place with perfect sequence
# what about linker A?
linkerA = lapply(fq$Sequence, substr, 57, 78)
100*table(unlist(lapply(linkerA, function(x) sum(strsplit(x, '')[[1]] != strsplit('ATCCACGTGCTTGAGACTGTGG', '')[[1]]))))/length(fq$Sequence)
#           0           1           2           3           4           5           6           7           8           9 
# 68.52583604  2.14748375  0.67971846  0.79197584  1.12257384  2.18172225  1.38862384  0.64660253  0.34182373  0.39795243 
#          10          11          12          13          14          15          16          17          18          19 
#  0.34743660  0.49337120  0.80993702  1.22023776  1.60752574  3.72919029  3.13534873  2.96527879  3.15387120  3.13029715 
#          20          21          22 
#  0.75493091  0.33003671  0.09822521
# 68% have linker A in the right place with perfect sequence
# let's investigate slippage
linkerAp1 = lapply(fq$Sequence, substr, 58, 79)
100*table(unlist(lapply(linkerAp1, function(x) sum(strsplit(x, '')[[1]] != strsplit('ATCCACGTGCTTGAGACTGTGG', '')[[1]]))))/length(fq$Sequence)
# no evidence of slippage
linkerAm1 = lapply(fq$Sequence, substr, 56, 77)
100*table(unlist(lapply(linkerAm1, function(x) sum(strsplit(x, '')[[1]] != strsplit('ATCCACGTGCTTGAGACTGTGG', '')[[1]]))))/length(fq$Sequence)
# no evidence of slippage
bc1 = lapply(fq$Sequence, substr, 79, 86)
bc18 = bc1[which(nchar(bc1)==8)]
barc1s = table(unlist(bc18))
barc1s = barc1s[order(barc1s, decreasing =T)]
# CTTTGGTC TTATTCTG GTGCTAGC GGGCGATG GGGTAGCG TCTATTAC CTGTCCCG GGTGGAGC TGGTATAC GTCGCGCG TTCCGATC TGGCGCGC ATATTGGC ACGCCGGC 
#    10877     3962     3241     3092     3015     2876     2833     2581     2217     2124     1986     1912     1833     1831 
# GCTCGCGG TAATACGC GACCTTTC TGCTTGGG AATTTCTC ATAAGCTC CGTCTAGG CACAATTG GCAAATTC TGGGCATC TTACCTCG TCTTAATC ATTCATGG CGTGGTTG 
#     1817     1636     1520     1504     1424     1413     1311     1273     1155     1144     1107     1070     1052     1017 
# CTGAAAGG TCGTTTCG TAAATATC GGTTCTTC GTGGGTTC TTCGCTAC GTACGACT ACGGACTC AGCGAAAC GAGACTGT TGCGATCG CTGCTTTG GAATAATG GCTGCATG 
#     1003      947      699      638      622      600      580      563      550      543      540      532      526      471 
# TTTGGTCT GTCATATG TTCATCGC GACAAAGC TTTTTTTT GCTATGCG ACGTTAAC CCATCTTG 
#      470      463      413      402      380      352      313      309
# bc_data_v2
# steep drop off after 50 but should only be 48... ok what about poly T, brings down to 49
bc2 = lapply(fq$Sequence, substr, 49, 56)
bc28 = bc2[which(nchar(bc2)==8)]
barc2s = table(unlist(bc28))
barc2s = barc2s[order(barc2s, decreasing =T)]
# bc_data_v1
bc3 = lapply(fq$Sequence, substr, 11, 18)
bc38 = bc3[which(nchar(bc3)==8)]
barc3s = table(unlist(bc38))
barc3s = barc1s[order(barc3s, decreasing =T)]
# bc_data_v1
```

# Basic QC on processed reads 

```{.r}
library('data.table')
WD = '/home/patrick/test/08_feature_count/'
expr = fread(paste0(WD, 'watch-sublibrary-1_S1.deduplicated.bam.featureCounts.tsv'))
dim(expr)
cSums = apply(expr[,-1], 2, sum)
sort(table(cSums))
