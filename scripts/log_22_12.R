# do var calling: barcode_headAligned_anno.bam - has all you need
#
# for expression matrix, filtered has only the nuclei with enough barcodes, while unfiltered has all the barcodes
library('Matrix')
library('future')
library('future.apply')
library('ineq')
library('ggplot2')
library('RColorBrewer')
library('ggplot2')
library('data.table')
plan(multicore, workers=20)
setDTthreads(20L)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')
setwd('~/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered/pipeline_output_data')
# reading in the data
# {{{
expr = readMM('DGE.mtx')
exprDf = data.table(base::as.matrix(t(expr)))
genes = read.csv('all_genes.csv')
meta = read.csv('cell_metadata.csv')
cellN = paste(meta$sample, meta$bc1_well, meta$bc2_well, meta$bc3_well, sep = '_')
colnames(exprDf) = cellN
rownames(exprDf) = genes$gene_id
exprDfT = exprDf[, lapply(.SD, mean), keyby = list(genes$gene_name), .SDcols = colnames(exprDf)[1:5000]]
# }}}

# gini calculation and plotting
# {{{
# calculate gini inequality for all cells
giniV = data.frame('Total_reads' = future_apply(exprDf, 2, sum, future.chunk.size = 1000L), 'Gini_inequality' = future_apply(exprDf, 2, Gini, future.chunk.size = 1000L))
giniV = data.frame(giniV, sample = meta$sample)
giniV$sample = gsub('s118', 's118\n38,451,432 reads', giniV$sample)
giniV$sample = gsub('s44', 's44\n7,181,054 reads', giniV$sample)
giniV$sample = gsub('s76', 's76\n16,240,929 reads', giniV$sample)
giniV$sample = gsub('s9', 's9\n13,536,611reads', giniV$sample)
giniV$sample = factor(giniV$sample, levels = c('s9\n13,536,611reads', 's44\n7,181,054 reads', 's76\n16,240,929 reads', 's118\n38,451,432 reads'))
# giniV = giniV[order(giniV$Total_reads, giniV$Gini_inequality, decreasing = c(TRUE, TRUE),]
pdf('all_gini.pdf')
print(ggplot(giniV, aes(x = Total_reads, y = Gini_inequality, color= sample))+
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_color_brewer(palette = 'Set1') +
    scale_x_continuous(limits = c(0, 200000), breaks = c(seq(0, 50000, 10000), seq(100000, 200000, 50000))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = 'Total number of unique reads per cell', y = 'Gini inequality', title = 'Effect of reads per cell on Gini inequality', subtitle = '15k reads is 833 cells, 20k reads is 567 cells, 25k reds is 409 cells')
)
dev.off()

pdf('all_gini_violin.pdf')
print(ggplot(giniV, aes(x = sample, y = Gini_inequality, fill = sample ))+
    geom_violin() +
    theme_bw() +
    scale_y_continuous(limits = c(0.99, 1), breaks = seq(.99, 1, .005)) +
    guides(fill="none") +
    labs(x='', y = 'Gini inequality', title = 'Violin plot of Gini inequality by sample', subtitle = 'Gini inequality > 0.99') +
    scale_fill_brewer(palette = 'Set1') 
)
dev.off()

giniVTrim = giniV[giniV$Gini_inequality<0.98,]
pdf('all_gini_Trim_violin.pdf')
print(ggplot(giniVTrim, aes(x = sample, y = Gini_inequality, fill = sample ))+
    geom_violin() +
    theme_bw() +
#    scale_y_continuous(limits = c(0.99, 1), breaks = seq(.99, 1, .005)) +
    guides(fill="none") +
    labs(x='', y = 'Gini inequality', title = 'Violin plot of Gini inequality by sample', subtitle = 'Gini inequality < 0.98') +
    scale_fill_brewer(palette = 'Set1') 
)
dev.off()
# }}}

# first filter will be to subselect to cells which have at least 5 reads in more than 100 genes
exprH = exprDf[, future_lapply(.SD, function(x) length(which(x>5)), future.chunk.size = 5000L)]
# plot of the number of cells which have more than 5 reads for X number of genes
# {{{
n = 400
plotG = data.frame('Number_of_genes' = seq(1,n), 'Number_of_cells' =unlist(lapply(seq(1,n), function(x) length(which(exprH>x)))))
pdf('all_gene_cutoff.pdf')
print(ggplot(plotG, aes(x = Number_of_genes, y = Number_of_cells))+
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_continuous(limits = c(0, 400), breaks = c(seq(0, 100, 20), seq(150, 400, 50))) +
    scale_y_continuous(limits = c(0, 35000), breaks = c(seq(10000, 30000, 10000), seq(0, 10000, 1000))) +
    labs(x = 'Number of genes with more than 5 reads required', y = 'Number of cells', title = 'Surviving cells after requiring 5 reads for X number of genes')
)
dev.off()
# }}}
#  MALAT1 MT-RNR1 MT-RNR2  KCNIP4  RBFOX1  SNHG14   OPCML   CSMD1    NRG3   NRXN3   LSAMP    DLG2   NRXN1    MEG3  DLGAP1 CNTNAP2   DPP10   PCDH9 
#  756023  481220  280447  230980  190549  141738  116305  116055   96587   96100   93909   91878   90950   89295   85756   85476   83597   82700 
#  ADGRB3  LRRTM4   ASIC2 FAM155A  ANKS1B    KAZN  DLGAP2   PTPRD    SYT1   FGF14   LRP1B   CADM2 
#   78580   75579   74041   72975   72456   70241   70231   69448   68638   67865   67858   67332 
#   sum(sort(z, decreasing = T))
# [1] 35409135
exprDfT = exprDf[,as.vector(exprH>200), with = F]

# second filter will be to genes that have more than 10 reads in at least 10% of cells
exprDfTa = transpose(exprDfT)
exprDfT = exprDfTa
exprI  = exprDfT[, future_lapply(.SD, function(x) length(which(x>10)), future.chunk.size = 3000L)]
# {{{
n = 841
plotG = data.frame('Number_of_cells' = seq(1,n), 'Number_of_genes' =unlist(lapply(seq(1,n), function(x) length(which(exprI>x)))))
library('ggplot2')
pdf('all_cell_cutoff.pdf')
print(ggplot(plotG, aes(x = Number_of_cells, y = Number_of_genes))+
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_continuous(limits = c(0, 850), breaks = c(seq(0, 200, 20), seq(250, 850, 50))) +
    scale_y_continuous(limits = c(0, 6100), breaks = c(seq(0, 2000, 250), seq(2500, 6000, 1000))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(y = 'Number of genes', x = 'Number of cells', title = 'Surviving genes after requiring 10 reads for X number cells')
)
dev.off()
# }}}

exprDfTa = transpose(exprDfT)
exprDfT = exprDfTa
exprDfTa = exprDfT[as.vector(exprI>167)]
colnames(exprDfTa) = cellN[as.vector(exprH>200)]
rownames(exprDfTa) = genes$gene_name[as.vector(exprI>167)]

exprDfTad = as.data.frame(exprDfTa)
colnames(exprDfTad) = cellN[as.vector(exprH>200)]
rownames(exprDfTad) = genes$gene_name[as.vector(exprI>167)]

rownames(exprDfwt) =exprDfwt$Gene
corT = cor(t(data.frame(exprDfTad)), method='spearman')
pdf('test.pdf')
heatmap(corT)
dev.off()

sanI = apply(exprDf, 1, mean)
exprDfTS = transpose(exprDfT)
exprDfT = exprDfTS

out = future_lapply(seq(1000, 5000, 1000), function(x){
    exprDfTT = exprDf[, 1:x]
    exprDfTTs = transpose(exprDfTT)
    exprDfTT = exprDfTTs
    exprS  = exprDfTT[, future_lapply(.SD, mean, future.chunk.size = 5000L)]
    return(length(which(unlist(exprS)>1)))
})

sanityP = data.frame(Cell_cutoff = seq(500,6000,500), Gene_number = unlist(out))
pdf('sanity_cutoff.pdf')
print(ggplot(sanityP, aes(x = Cell_cutoff, y = Gene_number))+
    geom_line() +
    geom_point() +
    theme_bw() +
    geom_vline(xintercept = 4000, color='red') +
    scale_x_continuous(limits = c(0, 6000), breaks = seq(0, 6000, 500)) +
    scale_y_continuous(limits = c(1500, 7500), breaks = seq(1500, 7500, 500)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(y = 'Number of genes', x = 'Number of cells', title = 'Surviving genes after requiring on average of 1 read across X number of cells')
)
dev.off()

# going to move forward with top 4000 cells and 2527 genes
exprDfC = exprDf[,1:4000]
exprDfCs = transpose(exprDfC)
exprDfC = exprDfCs
exprS  = exprDfC[, future_lapply(.SD, mean, future.chunk.size = 5000L)]
san = which(unlist(exprS)>1)
exprDfCs = transpose(exprDfC)
exprDfC = exprDfCs
exprDfC = exprDfC[san]
rownames(exprDfC) = genes$gene_name[san]
colnames(exprDfC) = cellN[1:4000]


# next, need to normalize data by library size
exprDfCN = log2(transpose(transpose(exprDfC) / unlist(exprDfC[, lapply(.SD, mean)]))+1)
rownames(exprDfCN) = genes$gene_name[san]
colnames(exprDfCN) = cellN[1:4000]
fwrite(exprDfCN, file = 'norm_log_4000nuc_2527gene_mat.csv', sep = ',')
# stop 
exprDfCN = fread('norm_log_4000nuc_2527gene_mat.csv')
# now hierarchical clustering
library('flashClust')
library('parallelDist')
library('cluster')
exprD = parallelDist(as.matrix(exprDfCN), method = 'euclidean', threads = 15)
exprC = flashClust(exprD, method = 'ward')
pdf('clust.pdf')
plot(exprC)
dev.off()
nClust=seq(2,50,1)
clusts = lapply(nClust, function(k) {
    return(cutree(exprC, k=k))
})

# 1) elbow plot
# function to calculate within sum of squared
calc_SS <- function(df) {
    sum(as.matrix(dist(df)^2)) / (2 * nrow(df))
}
# total cluster within sum of squared
plan(multicore, workers=20)
tclustwss= future_lapply(clusts, function(clust) {
    temp=c()
    for(ea in unique(clust)){
        ind = which(clust==ea)
        temp = c(temp, calc_SS(as.matrix(exprD)[ind, ind]))
    }
    return(sum(temp))
})
tclustwss = data.frame(nClust, unlist(tclustwss))

# 2) sillhoutte plot
sil = future_lapply(nClust, function(k) {
    distClust = cutree(exprC, k=k)
    return(silhouette(x = distClust, dmatrix = as.matrix(exprD)))
})
silW = future_lapply(sil, function(ea) {
    mean(as.matrix(ea)[,3])
})
silW = data.frame(nClust, unlist(silW))

# 3) gap statistic plot
hclusCut <- function(x, k, d.meth = "ward", ...) {
   list(cluster = cutree(hclust(exprD, d.meth, ...), k=k))
}
gapS = clusGap(as.matrix(exprD), FUN = hclusCut, K.max = 40, B = 10)

# ploting the results
pdf('hClust_k_choice.pdf')
plot(tclustwss, 
    xlab = 'Number of clusters', 
    ylab = 'Total within-clusters sum of squares', 
    main = 'Clustering choice dictated by elbow plot second derivative maxima',
    sub = 'Elbow at k=7,14',
    type = 'o')
plot(silW, 
    xlab = 'Number of clusters', 
    ylab = 'Average silhouette width', 
    main = 'Clustering choice dictated by silhouette width maxima',
    sub = 'Local maxima at k=12',
    type = 'o')
plot(tclustwss[1:50,], 
    xlab = 'Number of clusters', 
    ylab = 'Total within-clusters sum of squares', 
    main = 'Clustering choice dictated by elbow plot second derivative maxima',
    sub = 'Elbow at k=7,14',
    type = 'o')
plot(silW[1:50,], 
    xlab = 'Number of clusters', 
    ylab = 'Average silhouette width', 
    main = 'Clustering choice dictated by silhouette width maxima',
    sub = 'Local maxima at k=12',
    type = 'o')
plot(sil[[1]])
plot(sil[[2]])
plot(sil[[3]])
plot(sil[[4]])
plot(sil[[5]])
plot(sil[[6]])
plot(sil[[7]])
plot(sil[[8]])
plot(sil[[9]])
plot(sil[[10]])
dev.off()

# now with spearman's cor as dist
exprD = as.dist(1-cor(exprDfCN, method = 'spearman'))
exprC = flashClust(exprD, method = 'ward')
pdf('clust_spearman_cor.pdf')
plot(exprC)
dev.off()
nClust=seq(2,50,1)
clusts = lapply(nClust, function(k) {
    return(cutree(exprC, k=k))
})

# 1) elbow plot
# function to calculate within sum of squared
calc_SS <- function(df) {
    sum(as.matrix(dist(df)^2)) / (2 * nrow(df))
}
# total cluster within sum of squared
plan(multicore, workers=20)
tclustwss= future_lapply(clusts, function(clust) {
    temp=c()
    for(ea in unique(clust)){
        ind = which(clust==ea)
        temp = c(temp, calc_SS(as.matrix(exprD)[ind, ind]))
    }
    return(sum(temp))
})
tclustwss = data.frame(nClust, unlist(tclustwss))

# 2) sillhoutte plot
sil = future_lapply(nClust, function(k) {
    distClust = cutree(exprC, k=k)
    return(silhouette(x = distClust, dmatrix = as.matrix(exprD)))
})
silW = future_lapply(sil, function(ea) {
    mean(as.matrix(ea)[,3])
})
silW = data.frame(nClust, unlist(silW))

# ploting the results
pdf('hClust_k_choice_spearman.pdf')
plot(tclustwss, 
    xlab = 'Number of clusters', 
    ylab = 'Total within-clusters sum of squares', 
    main = 'Clustering choice dictated by elbow plot second derivative maxima',
    sub = 'Elbow at k=7,14',
    type = 'o')
plot(silW, 
    xlab = 'Number of clusters', 
    ylab = 'Average silhouette width', 
    main = 'Clustering choice dictated by silhouette width maxima',
    sub = 'Local maxima at k=12',
    type = 'o')
plot(tclustwss[1:50,], 
    xlab = 'Number of clusters', 
    ylab = 'Total within-clusters sum of squares', 
    main = 'Clustering choice dictated by elbow plot second derivative maxima',
    sub = 'Elbow at k=7,14',
    type = 'o')
plot(silW[1:50,], 
    xlab = 'Number of clusters', 
    ylab = 'Average silhouette width', 
    main = 'Clustering choice dictated by silhouette width maxima',
    sub = 'Local maxima at k=12',
    type = 'o')
for(i in seq(1, 20)){
    plot(sil[[i]])
}
dev.off()

# going forward with spearman cor, k = 16
library('parallelDist')
library('flashClust')
library('ComplexHeatmap')
library('dendextend')
source("/home/patrick/code/git/GSEA_generic/GSEAfxsV3.r")
source('~/code/git/GSEA_generic/enrichment_functions.R')

exprDfCL = (2^exprDfCN)-1
cZero = apply(exprDfCL, 2, function(x) length(which(x==0)))
rownames(exprDfCL) = rownames(exprDfCN)
csum = colSums(exprDfCL)
rsum = rowSums(exprDfCL)
exprDfCL = data.frame(exprDfCL)[order(rsum, decreasing = TRUE),order(cZero, decreasing = FALSE)]
rownames(exprDfCL) = rownames(exprDfCN)
exprDfCLFull = exprDfCL

for(kInd in c(16)){
for(iGene in c(seq(1000, 2000, 500), 2527)){
    for(jCell in seq(2000, 4000, 500)){
        # subselection and getting distance and dendrogram
        # {{{
        exprDfCL = exprDfCLFull[1:iGene, 1:jCell]
        exprD = as.dist(1-cor(exprDfCL, method = 'spearman'))
        # exprD = parallelDist(t(as.matrix(exprDfCL)), method = 'euclidean')
        exprC = flashClust(exprD, method = 'ward')
        dendP = as.dendrogram(exprC)
        dendP = color_branches(dendP, 
                                k = kInd)
        #                        col = dendCol)
        clustO = cutree(exprC, k=kInd)
        indOrd = unique(clustO[match(labels(exprC), names(clustO))])
        # }}}
        # p-value determination
        # {{{
        out = future_lapply(indOrd, function(clu){
            ind = (clu == clustO)
            future_apply(exprDfCL, 1, function(x){
                wilcox.test(x[ind], x[!ind], alternative='greater')$p.value
            })
        })
        # }}}
        # getting top 5 differentially expressed genes for each clone for heatmap
        # {{{
        exprN = future_lapply(seq_along(indOrd), function(i){
            ind = order(out[[i]], decreasing = F)
            try = data.frame(rownames(exprDfCL)[ind], exprDfCL[ind,])
            try[1:5,]
        })
        exprN = do.call(rbind, exprN)
        rownames(exprN) = make.unique(exprN[,1])
        plot = exprN[,-1]
        plot = t(apply(plot, 1, scale))
        colnames(plot) = colnames(exprN)[-1]
        plot = t(plot)
        # }}}
        # plotting heatmap
        # {{{
        setwd('~/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered')
        pdf(paste0(iGene, 'genes_', jCell, 'cells_',kInd, 'k_', 'sn_heatmap.pdf'), height=7, width = 20)
        print(Heatmap(plot,
            name = 'Scaled log2\nexpression',
            cluster_rows = dendP,
            col = colorRamp2(breaks = c(-4, 0, 4),
                                     colors = c('#2166ac', '#f7f7f7', '#b2182b')
                                    ),
            row_split = length(unique(clustO)),
            row_dend_width = unit(2, 'cm'),
        #     row_title = c('Clone 1', 'Clone 5', 'Astrocyte', 'Clone 3', 'Unknown', 'Clone 4 : 2', 'Clone 4 : 1','Oligodendrocyte', 'Microglia', 'Neuron')[seq(10,1)],
        #     row_title_gp = gpar(col = dendCol, fontfamily = 'NimbusSan', cex = 1.5),
            row_title_rot = 0,
            column_title = NULL,
            cluster_columns = FALSE,
            column_names_rot = 45, 
            column_names_side = "top",
        #     column_names_gp = gpar(col = rep(dendCol, ea =5), fontfamily = 'NimbusSan', cex = 0.8),
            column_split = rep(seq(1,length(unique(clustO))), ea =5)
        #     right_annotation = row_ha,
        #     left_annotation = row_ha_malig
        ))
        dev.off()
        # }}}
        # now need to do do enrichments of highly enriched genes...
        # {{{
        tout = lapply(indOrd, function(i) {
            temp = sort((out[[i]]*21000), decreasing = FALSE)
            temp = names(temp)[temp<0.01]
            temp
        })
        names(tout) = paste0('cluster', seq_along(indOrd))
        enrichs = enrichment_man(tout, rownames(exprDfCL), '/home/shared/genesets/genesets_slim')
        rm('mySets')
        setwd('~/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered')
        write.csv(enrichs, file = paste0(iGene, 'genes_', jCell, 'cells_',kInd,'k_','sn_enrich.csv'), row.names = F)
        if(!(exists('broadSets'))){
            broadSets=getBroadSets("/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml")
        }
        enrichs = enrichment_man_broad(tout, rownames(exprDfCL), broadSets)
        setwd('~/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered')
        write.csv(enrichs, file = paste0(iGene, 'genes_', jCell, 'cells_',kInd, 'k_', 'sn_enrich_broad.csv'), row.names = F)
        # }}}
    }
}}

enrichments = read.csv('sn_enrich.csv')
enrichments_broad = read.csv('sn_enrich_broad.csv')
colnames(enrichments_broad) = colnames(enrichments)
enrichments = rbind(enrichments, enrichments_broad)
keepSet = c('MOSET6944', 'MOSET6955', 'MOSET6945',
            'MOSET7051', 'MOSET133', 'MOSET7',
            'MOSET8', 'MOSET6936', 'MOSET7026',
            'MOSET9', 'MOSET7016', 'MOSET6941',
            'MOSET7016', 'MOSET9', 'MOSET40',
            'MOSET9', 'MOSET148', 'MOSET7057',
            'MOSET9', 'MOSET126', 'MOSET6941',
            'MOSET6933', 'MOSET6939', 'MOSET39',
            'MOSET6', 'MOSET148', 'MOSET26', 
            'MOSET7051', 'MOSET6', 'MOSET7025')

# nameSet = c('Verhaak: mesenchymal subtype', 'Tesar: OPC', 'Turcan: up in IDH1 mut',
#             'Suva: astrocytoma program', 'Kelley: astrocyte', 'Phillips: up in proneural subtype',
#             'Bachoo: astrocytes', 'Foster: mitochondrial proteins', 'Meissner: CpG promoters\nH3K4me3/H3K27me3',
#             'Zeisel: ependymal', 'Kelley: ependymal', 'Kelley: choroid',
#             'Verhaak: proneural subtype', 'Noushmehr: up in proneural', 'HPRD: androgen-\nreceptor pathway',
#             'Barres: astrocyte', 'Zeisel: astrocyte', 'Noushmehr: down in proneural',
#             'Kelley: oligodendrocyte', 'Zeisel: oligodendrocyte', 'Barres: oligodendrocyte',
#             'LaManno: microglia', 'Kelley: microglia', 'Suva: microglia',
#             'Barres: neuron', 'Kelley: neuron', 'ABA: neuron',
#             'Butler: endothelial', 'Kelley: endothelial', 'LaManno: endothelial'
#             )
 enrichments = enrichments[match(keepSet,enrichments$SetID),]
nameSet = enrichments$SetName
enrichments = as.matrix(enrichments[,8:ncol(enrichments)])
enrichments[enrichments == 0] = 1E-62
enrichments = -log10(enrichments)
rownames(enrichments) = nameSet

# do cor distance to get all malig clones together
pdf('single_nuc_enrichments.pdf', height = 13, width = 10)
Heatmap(enrichments,
    name = '-log10 p-values',
    cluster_columns = TRUE,
    clustering_method_columns = 'ward.D',
    clustering_distance_columns = 'pearson',
    cluster_rows = TRUE,
    clustering_method_rows = 'ward.D',
    clustering_distance_rows = 'pearson',
    column_names_side = "top",
    column_names_centered = TRUE,
    column_names_gp = gpar(fontfamily = 'NimbusSan', cex = 0.8),
    row_title_gp = gpar(fontfamily = 'NimbusSan', cex = 1.5),
    row_dend_width = unit(2, 'cm'),
    row_title = NULL,
    column_names_rot = 45, 
    row_split = rep(seq(1, 10), ea=3),
#     left_annotation = rha,
#     right_annotation = row_ha,
    col = colorRamp2(breaks = c(0, 15, 66),
                            colors = c('#ffffff', '#fc9272', '#de2d26')
                            )
)
dev.off()

# do find modules
library("future")
library('data.table')
exprDfTG = exprDfT$genes
exprDfT = exprDfT[, -1, with =F]
exprH = exprDfT[, future_lapply(.SD, function(x) length(which(x>5)), future.chunk.size = 5000L)]
exprDfFT = exprDfT[, order(exprH, decreasing = TRUE)[1:5000], with = F]
exprDfT = exprDfFT[, lapply(.SD, sum), keyby = list(exprDfTG), .SDcols = colnames(exprDfFT)[1:5000]]
exprDfFM = transpose(exprDfT[,-1, with=F])
nonZ = exprDfFM[,future_lapply(.SD, function(x) length(which(x>0)), future.chunk.size = ncol(exprDfFM)/5)]
exprDfFM = exprDfFM[, nonZ>50, with = F]
exprDfFM = transpose(exprDfFM)
rownames(exprDfFM) = exprDfT$genes[nonZ>50]
colnames(exprDfFM) = colnames(exprDfFT)
exprDfFMF = t(t(data.frame(exprDfFM)) / unlist(exprDfFM[, lapply(.SD, sum)]))
exprOut = data.frame(Gene = rownames(exprDfFM), exprDfFMF)
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
exprOut = data.frame(fread('/mnt/bdata/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered/normalized_matrices/19127_genes_5000cells_orderedbynonzero.csv'))
exprOut[,-1] =exprOut[,-1] *1E4
setwd('/mnt/bdata/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered/')
FindModules(
    projectname='19127_genes_5000cells',
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
    signumvec=c(seq(0.5, 0.9, 0.1), 0.95, 0.99),
    minsizevec=c(15, 10, 8, 5,3),
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

# new
# do copykat call, try, look for 1p19q
# need code that was used to run copykat and unadjusted count data
# import raw read counts
library('Matrix')
library('future')
library('future.apply')
library('ineq')
library('ggplot2')
library('RColorBrewer')
library('ggplot2')
library('data.table')
plan(multicore, workers=20)
setDTthreads(20L)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')
setwd('~/@patrick/parsebio/dec_run/03_parse_pipe_proc/all-well/DGE_unfiltered')
expr = readMM('DGE.mtx')
exprDf = transpose(data.table(as.matrix(expr)))
genes = read.csv('all_genes.csv')
meta = read.csv('cell_metadata.csv')
cellN = paste(meta$sample, meta$bc1_well, meta$bc2_well, meta$bc3_well, sep = '_')
colnames(exprDf) = cellN
rownames(exprDf) = genes$gene_id



# run copykat
library('copykat')
library('data.table')
library('ggplot2')
library('future')
library('future.apply')
library('qvalue')
library('dendextend')
library('RColorBrewer')
library('grid')
library('gridExtra')
library('limma')
library('ComplexHeatmap')
plan(multicore, workers=20)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')

# will include top 4000 cells by read count. 
# {{{
# questions: 
# what does low and high dr mean? cutoffs to only keep genes that are expressed in the low.dr to up.dr fraction of cells
# what is ngene.chr? - number of genes necessary to represent a chromosome
# what is win.size? - window size on which to make call, in mb
# distance - pearson/spearman better for noisy data, euclidean better for high read data
# what is Ks.cut - sensititivity of algorithm, default is 0.1, higher decreases sensitivity
# output.seg - output segmentation file for IGV visualization
# plot.genes - output plot of copy number calls
source('/opt/copykat/R/copykat.R')
exprDfT = as.data.frame(temp)
rownames(exprDfT) = exprDfT[,1]
exprDfT = exprDfT[,-1]
copyDat = copykat(rawmat = exprDfT[, 1:2000], 
                    id.type = 'Symbol',
                    cell.line = 'no',
                    ngene.chr = 20, 
                    LOW.DR = 0.05,
                    UP.DR = 0.3,
                    win.size = 10, 
                    norm.cell.names = NULL,
                    KS.cut = 0.2, 
                    sam.name = 'copykat_parse', 
                    distance = 'spearman',
                    output.seg = FALSE,
                    plot.genes = TRUE,
                    genome = 'hg20',
                    n.cores = 10)

copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", genome="hg20",n.cores=1)
# }}}
