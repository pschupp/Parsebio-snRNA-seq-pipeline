library('Matrix')
library('future')
library('future.apply')
library('ineq')
library('ggplot2')
library('RColorBrewer')
library('ggplot2')
library('data.table')
library("future")
library('data.table')
plan(multicore, workers=20)
setDTthreads(20L)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')
setwd('~/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered/pipeline_output_data')
# reading in the data
expr = readMM('DGE.mtx')
exprDf = data.table(base::as.matrix(t(expr)))
genes = read.csv('all_genes.csv')
meta = read.csv('cell_metadata.csv')
cellN = paste(meta$sample, meta$bc1_well, meta$bc2_well, meta$bc3_well, sep = '_')
colnames(exprDf) = cellN
rownames(exprDf) = genes$gene_id

exprDfT = exprDf[, lapply(.SD, mean), keyby = list(genes$gene_name), .SDcols = colnames(exprDf)[1:5000]]
exprDfTG = exprDfT$genes
exprDfT = exprDfT[, -1, with =F]
exprH = exprDfT[, future_lapply(.SD, function(x) length(which(x>5)), future.chunk.size = 5000L)]
exprDfFT = exprDfT[, order(exprH, decreasing = TRUE)[1:5000], with = F]
exprDfT = exprDfFT[, lapply(.SD, sum), keyby = list(exprDfTG), .SDcols = colnames(exprDfFT)[1:5000]]
exprDfFM = transpose(exprDfT[,-1, with=F])
nonZ = exprDfFM[,future_lapply(.SD, function(x) length(which(x>0)), future.chunk.size = ncol(exprDfFM)/5)]
exprDfFM = exprDfFM[, nonZ>500, with = F]
exprDfFM = transpose(exprDfFM)
rownames(exprDfFM) = exprDfT$genes[nonZ>500]
colnames(exprDfFM) = colnames(exprDfFT)
exprDfFMF = t(t(data.frame(exprDfFM)) / unlist(exprDfFM[, lapply(.SD, sum)]))
exprOut = data.frame(Gene = exprDfTG[nonZ>50], exprDfFMF)
exprOut = exprOut[, c(1, order(gsub('_.*', '', colnames(exprOut)[-1]))+1)]
fwrite(exprOut, file = '/mnt/bdata/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered/normalized_matrices/19127_genes_5000cells_orderedbynonzero.csv', row.names=FALSE, sep = ',')
#### stop
library("WGCNA")
library("svMisc")
library("Biobase")
library("lattice") ## 3D plots
library("qvalue")
library("ellipse") ## plotcorr function
library("purrr") ## required for merging by CC
library("HiClimR") # for fastCor
library("RColorBrewer") 
library("future")
library('data.table')
source('~/code/git/FindModules/FindModules.lint.par.R')
plan(multicore, workers=5)
setDTthreads(20L)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')
exprOut = data.frame(fread('/mnt/bdata/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered/normalized_matrices/9930_genes_5000cells_orderedbynonzero.csv'))
exprOut[,-1] =exprOut[,-1] *1E4
source('~/code/git/FindModules/FindModules.R')
setwd('/mnt/bdata/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered/')
FindModules(
    projectname='9930_genes_5000cells',
    expr=exprOut,
    geneinfo=c(1),
    sampleindex=seq(2,ncol(exprOut)),
    samplegroups=as.factor(gsub('_.*', '', colnames(exprOut))),
    subset=NULL,
    simMat=NULL,
    saveSimMat=FALSE,
    simType="Bicor",
    beta=1,
    overlapType="None",
    TOtype="signed",
    TOdenom="min",
    MIestimator="mi.mm",
    MIdisc="equalfreq",
    signumType="rel",
    iterate=TRUE,
    signumvec=c(seq(0.8, 0.9, 0.1), 0.95, 0.99,.99),
    minsizevec=c(25, 20, 15, 10, 8),
    signum=NULL,
    minSize=NULL,
    merge.by="ME",
    merge.param=0.8,
    export.merge.comp=T,
    ZNCcut=2,
    calcSW=FALSE,
    loadTree=FALSE,
    writeKME=TRUE,
    calcBigModStat=FALSE,
    writeModSnap=TRUE
)


# gene set enrichment
library('flashClust')
library('parallel')
source("~/code/git/GSEA_generic/GSEAfxsV3.r")
WD="/mnt/bdata/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered/19127_genes_5000cells_Modules"
setwd(WD)
print("Set working directory")
## To run enrichment analysis for our gene sets in all networks:
MyGSHGloop(kmecut1="topmodposbc",exclude="none",pvalcut1=NULL)
setwd(WD)
MyGSHGloop(kmecut1="topmodposfdr",exclude="none",pvalcut1=NULL)
setwd(WD)

## Read in Broad gene sets (MolSigDBv3): 
#print("Reading in 'broadSets'")
#
broadSets=getBroadSets("/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml")
print("Success reading in 'broadSets'")
## To run enrichment analysis for broad gene sets in all networks:
### Note that kmecut1 can equal "seed", "topmodposbc" (recommended), or "topmodposfdr".
print("Beginning loop BC")
BroadGSHGloop(kmecut1="topmodposbc",pvalcut1=NULL)
setwd(WD)
print("Beginning loop FD")
BroadGSHGloop(kmecut1="topmodposfdr",pvalcut1=NULL)
setwd(WD)

