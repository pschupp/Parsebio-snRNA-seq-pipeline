# prep environment
# {{{
library('data.table')
library('limma')
library('future')
library('future.apply')
library('RColorBrewer')
library('grid')
library('gridExtra')
library('ComplexHeatmap')
library('parallelDist')
plan(multicore, workers=20)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')
# }}}
# 1. get all mutation names for each clone
# {{{
loci=fread('~/@patrick/SF10711/clonality/pyclone/outputs/tables/loci.tsv')
mutMember = aggregate(loci$cluster_id, by = list(loci$mutation_id), FUN = max)
colnames(mutMember) = c('mutation', 'clone')
mutMember = mutMember[order(mutMember$clone),]
# }}}
# 2. get expression percentile of each gene, from both single-nucleus and bulk
# {{{
# bulk rna expression matrix, full normalized
# {{{
rna=fread('/mnt/bdata/@patrick/SF10711/rna.seq/expression_matrices/Normalized_read_counts_using_RUVg_ERCC_K10Factors.csv', drop=seq(2,6))
colnames(rna)[1]='Gene'
meanExpr=apply(rna[,-seq(1,6)], 1, var)
rna=rna[order(meanExpr, decreasing=T),]
rnarSums=apply(rna[,-1], 1, sum)
# rna=rna[-which(rnarSums<quantile(rnarSums,0.15)),]
rna=rna[-grep('^ERCC-', rna$Gene),]
rnarSumsN=apply(rna[,-1], 1, sum)
rna=as.data.frame(rna)
# rna=rna[,c(1, which(colnames(rna) %in% barOut$ampseq.sample))]
rnaScale=t(future_apply(rna[,-c(1)], 1, log2))
# correct line below and make sure that the sections line up as expected
colnames(rnaScale)=gsub('SF10711)', '', colnames(rnaScale))
rnaScale=data.frame(Gene=alias2SymbolTable(rna$Gene, species='Hs'), rnaScale)
rnaScale=data.frame(Gene=rna$Gene, rnaScale)
# }}}
# single nucleus expression matrix, full normalized
# {{{
snExpr = fread('~/@patrick/SF10711/sn.rna.seq/final_data/expr.class.name.ruvg.k10.unround.cpm.log2.combat.csv', sep =',', header=T, data.table =F)
snExprGene = alias2SymbolTable(snExpr$Gene, species = 'Hs')
# meanExpr = apply(snExpr[,-1], 1, mean)
# snExpr = snExpr[order(meanExpr, decreasing =T)[-c(1:6000)],]
snExpr[,-c(1,2)] = log2(snExpr[ ,-c(1,2)]+1)
inds = list(
    s17 = grep('N701|N702|N709', colnames(snExpr)),
    s53 = grep('N701|N702|N709', colnames(snExpr)),
    s93 = grep('N701|N702|N709', colnames(snExpr)),
    s117 = grep('N701|N702|N709', colnames(snExpr))
    )
snExpr = data.frame(Gene = snExprGene, future_lapply(inds, function(x) apply(snExpr[,x], 1, mean)))
snExpr = data.frame(do.call(cbind, snExpr))
# }}}
# limit bulk to certain sections
# {{{
secs = list(17,53, 93, 117)
names(secs) = paste0('s', c(17,53, 93, 117))
# the following determines how many variables surrounding the single nucleus section to include
surSamples = 2
inds = lapply(secs, function(x) order(abs(x -  as.numeric(gsub('.*_', '', colnames(rnaScale)))))[1:surSamples])
bulkExpr = data.frame(Gene = rnaScale$Gene, future_lapply(inds, function(x) apply(rnaScale[,x], 1, mean)))
# jointGenes = intersect(snExpr$Gene, rnaScale$Gene)
# snExpr = snExpr[match(jointGenes, snExpr$Gene),-1]
# bulkExpr = (1-(rank(apply(bulkExpr[match(jointGenes, bulkExpr$Gene),-1], 1, mean))/length(jointGenes)))*100
snExpr = (1-(rank(apply(apply(snExpr, 2, as.numeric), 1, mean))/nrow(snExpr)))*100
bulkExpr = (1-(rank(apply(bulkExpr[,-1], 1, mean))/nrow(bulkExpr)))*100
# exprPlot = data.frame(sn = snExpr, bulk = bulkExpr)
# }}}
# out = data.frame(gene = jointGenes, exprPlot)
mutGenes = gsub('_.*', '', mutMember$mutation)
outMut = data.frame(snExpr = snExpr[match(mutGenes, snExprGene)], bulkExpr = bulkExpr[match(mutGenes, rnaScale$Gene)])
mutMember = data.frame(mutMember, outMut)
# }}}
# 3. get distance from end and exon number, etc.
vep  = fread('~/@patrick/SF10711/figures/liftover_GRCh37_to_38/GRCh38_positions.tsv', skip = 85)
mutMember = data.frame(mutMember, chr = vep$V1[match(mutGenes, vep$V5)], start = vep$V2[match(mutGenes, vep$V5)])
gtf=fread('/home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf', skip=5)
colnames(gtf)=c("chr", "db", "type", "start", "end", "score", "strand", "frame", "extra")
gtf=gtf[which(gtf$type=="exon"),]
gtf.ENSG=gsub(".*gene_id\\s\"\\\"*|\\..*", "", gtf$extra)
gtf.type=gsub(".*gene_type\\s\"\\\"*|\".*", "", gtf$extra)
gtf.name=gsub(".*gene_name\\s\"\\\"*|\".*", "", gtf$extra)
gtfOut$gtf.name= alias2SymbolTable(gtfOut$gtf.name, species = 'Hs')
gtfOut$gtf.name[is.na(gtfOut$gtf.name)] = gtfOut$gtf.ENSG[is.na(gtfOut$gtf.name)]
# get genome in memory to extract exon sequences 
# {{{
library('Biostrings')
genome = readDNAStringSet(filepath = '/home/shared/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa', format = 'fasta')
genome@ranges@NAMES = gsub('\\s.*', '', genome@ranges@NAMES)
# test_2bit_out <- file.path(tempdir(), "test_out.2bit")
# rtracklayer::export.2bit(wheat_genome_seqs, test_2bit_out)
# }}}
exon_start_get = function(gene, loc){
    tempGtf = gtf[gtf.name == gene,]
    tempGtf = tempGtf[grep('Ensembl_canonical', tempGtf$extra),]
    if(nrow(tempGtf)==0){
        return(NA)
        next
    }
    if(tempGtf$strand[1] == '-'){
        incStrand = FALSE
    } else if(tempGtf$strand[1] == '+'){
        incStrand = TRUE
    } 
    tempGtf = tempGtf[(tempGtf$start <= loc & tempGtf$end >= loc),]
    if(nrow(tempGtf)==0){
        return(NA)
        next
    }
    if(incStrand){
        dist = loc - tempGtf$start
        strandDir = 'forward'
        seq = genome[[which(names(genome) == tempGtf$chr)]][tempGtf$start:tempGtf$end]
    } else if(!(incStrand)){
        dist = tempGtf$end - loc
        strandDir = 'reverse'
        seq = reverseComplement(genome[[which(names(genome) == tempGtf$chr)]][tempGtf$start:tempGtf$end])
    }
    return(c(dist, strandDir, tempGtf$chr, tempGtf$start, tempGtf$end, as.character(seq)))
}    
start_ex = lapply(seq_len(nrow(mutMember)), function(x) exon_start_get(mutGenes[x], mutMember$start[x]))
start_ex = do.call(rbind, start_ex)
mutMember = data.frame(mutMember, start_ex)
mutMember$clone = gsub('^0', 'clone 5', mutMember$clone)
mutMember$clone = gsub('^1', 'clone 2', mutMember$clone)
mutMember$clone = gsub('^2', 'clone 1', mutMember$clone)
mutMember$clone = gsub('^3', 'clone 4', mutMember$clone)
mutMember$clone = gsub('^4', 'clone 3', mutMember$clone)
mutMember = mutMember[, -5]
colnames(mutMember) = c('mutation', 'clone', 'snExpr', 'bulkExpr', 'SNP_location', 'dist_to_5prime', 'strandedness', 'exon_start', 'exon_end', 'sequence')
mutMember = mutMember[order(mutMember$clone),]
write.table(mutMember, file = 'SF10711_targets.csv', sep=',', row.names=F, quote=F)
