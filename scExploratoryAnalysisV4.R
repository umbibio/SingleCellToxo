library(tidyverse)
library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)


##### Markers from Ben Waldman
markers.TB <- read.xlsx('~/work/ToxoPlasmaGondiiR/Input/markers/Gene_modules_Waldman2019_KZ.xlsx', sheet = 5)
markers.G1SM <- read.xlsx('~/work/ToxoPlasmaGondiiR/Input/markers/Gene_modules_Waldman2019_KZ.xlsx', sheet = 4)

ID.Orthologs <- read.xlsx('~/work/ToxoPlasmaGondiiR/Input/ID_convert/convertIDs.xlsx')

markers.G1SM <- left_join(markers.G1SM, ID.Orthologs, by = c('GeneID' = 'TGME49ID')) %>% 
  dplyr::select(TGGT1ID, Marker) %>% dplyr::filter(!is.na(Marker)) %>% transmute(GeneID = TGGT1ID, Marker = Marker)

markers <- rbind(markers.G1SM, markers.TB) 


microneme.markers <- read.xlsx('~/work/ToxoPlasmaGondiiR/Input/markers/Microneme life cycle expression -ScRNAseq.xlsx')

abs.path <- "~/work/ToxoPlasmaGondiiR/Input/SingleCell/"
SW3.meta.rds <- readRDS(file = paste(abs.path, "SW3_meta_data.rds",  sep = ""))
SW3.sparse.rds <- readRDS(file = paste(abs.path, "SW3_sparse_expression.rds",  sep = ""))

M <- as.data.frame(as.matrix(SW3.sparse.rds))
Meta <- as.data.frame(as.matrix(SW3.meta.rds))
Meta <- Meta %>% mutate(Sample = rownames(Meta))

## Convert IDs to GT1
M.GT1 <- M[!is.na(match(rownames(M), ID.Orthologs$TGME49ID)),]
rownames(M.GT1) <- ID.Orthologs$TGGT1ID[match(rownames(M.GT1), ID.Orthologs$TGME49ID)]

Meta.72 <- Meta %>% dplyr::filter(Timepoint == 72)
M.72 <- M.GT1 %>% dplyr::select(colnames(M)[colnames(M.GT1) %in% Meta.72$Sample])

M.72.KO.D3 <- M.72 %>% dplyr::select(contains('KO.D3'))
M.72.KO.pH <- M.72 %>% dplyr::select(contains('KO'),-contains('D3'))
M.72.WT.D3 <- M.72 %>% dplyr::select(contains('WT.D3'))
M.72.WT.pH <- M.72 %>% dplyr::select(contains('WT'), -contains('D3'))
M.72.WT    <- M.72 %>% dplyr::select(contains('WT'))
M.72.D3    <- M.72 %>% dplyr::select(contains('D3'))

#### Do some Stats
## for each gene, identify # cells in which the gene is present
getAbundantGenes <- function(M){
  tmp <- M
  tmp[tmp == 0] <- FALSE
  tmp[tmp != 0] <- TRUE
  num.cells.expressing <- rowSums(tmp)
  ## for each gene, identify mean read number accross cells
  mean.reads <- rowSums(M) / rowSums(tmp)
  
  q1 <- quantile(num.cells.expressing)
  q2 <- quantile(mean.reads, na.rm = T)
  
  abundant.genes <- which(num.cells.expressing >= q1[3] & mean.reads >= q2[3])
  
  abundant.genes <- data.frame(GeneID = names(abundant.genes), 
                               num_cells_expressing = num.cells.expressing[abundant.genes],
                               ave_expr = mean.reads[abundant.genes])
  return(abundant.genes)
}

# abundant.genes.WT.D3 <- getAbundantGenes(M.72.WT.D3)
# abundant.genes.WT.pH <- getAbundantGenes(M.72.WT.pH)
# abundant.genes.WT <- getAbundantGenes(M.72.WT)
# abundant.genes <- getAbundantGenes(M.72)
# abundant.genes.D3 <- getAbundantGenes(M.72.D3)
# 
# write.xlsx(abundant.genes.WT.D3, '~/work/ToxoPlasmaGondiiR/Output/Tables/abundant_genes_WT_D3_more_stringent.xlsx')  
# write.xlsx(abundant.genes.WT.pH, '~/work/ToxoPlasmaGondiiR/Output/Tables/abundant_genes_WT_pH_more_stringent.xlsx') 
# write.xlsx(abundant.genes.WT, '~/work/ToxoPlasmaGondiiR/Output/Tables/abundant_genes_WT_more_stringent.xlsx') 
## for each cell, identify the number of genes expressed

getAbundantCells <- function(M, abundant.genes){
  tmp <- M
  tmp[tmp == 0] <- FALSE
  tmp[tmp != 0] <- TRUE
  
  num.genes.expressed <- colSums(tmp)
  
  ## for each cell, identify mean read number accross all genes
  mean.reads.cells <- colSums(M[abundant.genes, ]) / colSums(tmp[abundant.genes,])
  
  prop.of.abundant.genes <- colSums(tmp[abundant.genes,]) / length(abundant.genes)
  
  q3 <- quantile(num.genes.expressed)
  q4 <- quantile(mean.reads.cells, na.rm = T)
  
  ## take cells expression at least 20% of the abundant genes, expression high number of genes, with high coverage 
  abundant.cells <- which(num.genes.expressed >= q3[2] & mean.reads.cells >= q4[3] & 
                            prop.of.abundant.genes > 0.2)
  
  return(abundant.cells)
}

# abundant <- M.72.WT.D3[abundant.genes, abundant.cells]
# 
# tmp <- abundant
# tmp[tmp == 0] <- FALSE
# tmp[tmp != 0] <- TRUE
# 
# plot(colSums(abundant), colSums(tmp))
# S.O.abundant <- CreateSeuratObject(counts =abundant, min.cells = 3, min.features = 200)
# VlnPlot(S.O.abundant, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
# FeatureScatter(S.O.abundant, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


## Now take the intersection of abundant genes with cyclic genes
#cyclic.markers <- markers[markers$Marker %in% c('G1', 'SM'), ]
#cyclic.abundant.markers <- inner_join(abundant.genes.WT.D3, cyclic.markers, by = 'GeneID')
#abundant.cyclic <- M.72.WT.D3[cyclic.abundant.markers$GeneID, abundant.cells]
#abundant.cyclic <- M.72.WT.D3[cyclic.abundant.markers$GeneID, ] ## take all cells
#tmp <- abundant.cyclic
#tmp[tmp == 0] <- FALSE
#tmp[tmp != 0] <- TRUE

# abundant.markers.WT.pH <- inner_join(abundant.genes.WT.pH, markers, by = 'GeneID')
# abundant.markers.WT.pH <- M.72.WT.pH[unique(abundant.markers.WT.pH$GeneID), ] ## take all cells
# 
# abundant.markers.WT.D3 <- inner_join(abundant.genes.WT.D3, markers, by = 'GeneID')
# abundant.markers.WT.D3 <- M.72.WT.D3[unique(abundant.markers.WT.D3$GeneID), ] ## take all cells
# 
# abundant.markers.WT <- inner_join(abundant.genes.WT, markers, by = 'GeneID')
# abundant.markers.WT <- M.72.WT[unique(abundant.markers.WT$GeneID), ] ## take all cells
# 
# abundant.markers <- inner_join(abundant.genes, markers, by = 'GeneID')
# abundant.markers <- M.72[unique(abundant.markers$GeneID), ] ## take all cells
# 
# abundant.markers.D3 <- inner_join(abundant.genes.D3, markers, by = 'GeneID')
# abundant.markers.D3 <- M.72.D3[unique(abundant.markers.D3$GeneID), ] ## take all cells

# tmp <- abundant.markers
# tmp[tmp == 0] <- FALSE
# tmp[tmp != 0] <- TRUE
# 
# plot(colSums(abundant.markers), colSums(tmp))

S.O.WT.D3 <- CreateSeuratObject(counts = M.72.WT.D3, min.cells = 3, min.features = 200)
VlnPlot(S.O.WT.D3, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.WT.D3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.WT.D3 <- subset(S.O.WT.D3, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 )

#WT.pH
S.O.WT.pH <- CreateSeuratObject(counts = M.72.WT.pH, min.cells = 3, min.features = 200)
VlnPlot(S.O.WT.pH, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.WT.pH, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.WT.pH <- subset(S.O.WT.pH, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 )

##KO.D3
S.O.KO.D3 <- CreateSeuratObject(counts = M.72.KO.D3, min.cells = 3, min.features = 200)
VlnPlot(S.O.KO.D3, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.KO.D3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.KO.D3 <- subset(S.O.KO.D3, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 )

#KO.pH
S.O.KO.pH <- CreateSeuratObject(counts = M.72.KO.pH, min.cells = 3, min.features = 200)
VlnPlot(S.O.KO.pH, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.KO.pH, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.KO.pH <- subset(S.O.KO.pH, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 )

prep_S.O <- function(S.O){
  S.O <- NormalizeData(S.O, normalization.method = "LogNormalize", scale.factor = 10000)
  S.O <- FindVariableFeatures(S.O, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(S.O)
  S.O <- ScaleData(S.O, features = all.genes)
  S.O <- RunPCA(S.O, features = VariableFeatures(object = S.O))
  S.O <- FindNeighbors(S.O, dims = 1:13)
  S.O <- FindClusters(S.O, resolution = 0.1)
  S.O <- RunUMAP(S.O, dims = 1:13)
  return(S.O)
}

S.O.WT.D3 <- prep_S.O(S.O.WT.D3)
S.O.WT.pH <- prep_S.O(S.O.WT.pH)
S.O.KO.D3 <- prep_S.O(S.O.KO.D3)
S.O.KO.pH <- prep_S.O(S.O.KO.pH)
S.O.abundant <- prep_S.O(S.O.abundant)

S.O.abundant.markers.WT.D3 <- prep_S.O(S.O.abundant.markers.WT.D3)
S.O.abundant.markers.WT.pH <- prep_S.O(S.O.abundant.markers.WT.pH)
S.O.abundant.markers.WT <- prep_S.O(S.O.abundant.markers.WT)
S.O.abundant.markers <- prep_S.O(S.O.abundant.markers)
S.O.abundant.markers.D3 <- prep_S.O(S.O.abundant.markers.D3)

DimPlot(S.O.WT.D3, reduction = "umap")
DimPlot(S.O.KO.D3, reduction = "umap")
DimPlot(S.O.WT.pH, reduction = "umap")
DimPlot(S.O.KO.pH, reduction = "umap")
DimPlot(S.O.abundant, reduction = "umap")

DimPlot(S.O.abundant.markers.WT.D3, reduction = "umap")
DimPlot(S.O.abundant.markers.WT.pH, reduction = "umap")
DimPlot(S.O.abundant.markers.WT, reduction = "umap")
DimPlot(S.O.abundant.markers, reduction = "umap")
DimPlot(S.O.abundant.markers.D3, reduction = "umap")
#top10 <- head(VariableFeatures(S.O), 10)

#plot1 <- VariableFeaturePlot(S.O)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))

#print(S.O[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(S.O, dims = 1:2, reduction = "pca")
#DimPlot(S.O, reduction = "pca")
#DimHeatmap(S.O, dims = 1, cells = 500, balanced = TRUE)

#S.O <- JackStraw(S.O, num.replicate = 100)
#S.O <- ScoreJackStraw(S.O, dims = 1:20)
#JackStrawPlot(S.O, dims = 1:20)
#ElbowPlot(S.O)



## Get normalized expression data
getNormExpr <- function(S.O, markers){
  expr.norm <- as.matrix(S.O[["RNA"]]@data)
  rownames(expr.norm) <- gsub('-', '_',rownames(expr.norm))
  expr.norm <- as.data.frame(expr.norm) 
  expr.norm <- expr.norm %>% mutate(GeneID = rownames(expr.norm))
  expr.norm <- expr.norm %>% gather(key = Sample, value = expr, -GeneID)
  expr.norm <- expr.norm  %>% group_by(GeneID) %>% 
    mutate(quantile = ifelse(expr <= quantile(expr)[2], 1, 
                             ifelse(expr <= quantile(expr)[3], 2, 
                                    ifelse(expr <= quantile(expr)[4], 3,4))))
  
  expr.norm <- left_join(expr.norm, markers, by = 'GeneID')
  
  return(expr.norm)
}

## Get normalized expression Manually
getZscores <- function(M, markers){
  expr.norm <- M
  rownames(expr.norm) <- gsub('-', '_',rownames(expr.norm))
  expr.norm <- as.matrix(expr.norm) 
  r.m <- matrix(matrixStats::rowMeans2(expr.norm), nrow = nrow(expr.norm), ncol = ncol(expr.norm))
  r.sd <- matrix(matrixStats::rowSds(expr.norm), nrow = nrow(expr.norm), ncol = ncol(expr.norm))
  expr.norm <- (expr.norm - r.m) / r.sd
  expr.norm <- as.data.frame(expr.norm) 
  expr.norm <- expr.norm %>% mutate(GeneID = rownames(expr.norm))
  expr.norm <- expr.norm %>% gather(key = Sample, value = expr, -GeneID)
  expr.norm <- expr.norm  %>% group_by(GeneID) %>% 
    mutate(quantile = ifelse(expr <= quantile(expr)[2], 1, 
                             ifelse(expr <= quantile(expr)[3], 2, 
                                    ifelse(expr <= quantile(expr)[4], 3,4))))
  
  expr.norm <- left_join(expr.norm, markers, by = 'GeneID')
  
  return(expr.norm)
}

expr.norm.WT.D3      <- getNormExpr(S.O.WT.D3, markers)
expr.norm.WT.pH      <- getNormExpr(S.O.WT.pH, markers)
expr.norm.KO.D3      <- getNormExpr(S.O.KO.D3, markers)
expr.norm.KO.pH      <- getNormExpr(S.O.KO.pH, markers)
expr.norm.abundant   <- getNormExpr(S.O.abundant, markers)
expr.zscore.abundant <- getZscores(abundant, markers)


expr.norm.abundant.markers.WT   <- getNormExpr(S.O.abundant.markers.WT, markers)
expr.norm.abundant.markers.WT.D3   <- getNormExpr(S.O.abundant.markers.WT.D3, markers)
expr.norm.abundant.markers.WT.pH   <- getNormExpr(S.O.abundant.markers.WT.pH, markers)
expr.norm.abundant.markers   <- getNormExpr(S.O.abundant.markers, markers)
expr.norm.abundant.markers.D3   <- getNormExpr(S.O.abundant.markers.D3, markers)

getCellIdentity <- function(expr.norm, method = ''){
  cell.identity <- expr.norm %>% ungroup() %>% dplyr::select(Sample, expr, quantile, Marker) %>% distinct() 
  cell.identity <- cell.identity %>% dplyr::filter(!is.na(Marker))
  cell.identity <- cell.identity %>% group_by(Sample, Marker) %>% 
    summarise(mean.expr = mean(expr),
              sd.expr = sd(expr),
              max.quantile = max(quantile), 
              prop.expressed = sum(quantile >= 2) / n(),
              prop.highly.expressed = sum(quantile >= 3) / n(),
              prop.repressed = sum(quantile < 2) / n())
  #,zscore = sum(expr > 0) / n())
  ## Using proportion of highly expressed for identity
  if(method == 'prop.highly.expressed'){
    cell.identity <- cell.identity %>% dplyr::select(Sample, Marker, prop.highly.expressed) %>%
      spread(key = Marker, value = prop.highly.expressed) 
  }else if(method == 'prop.expressed'){
    cell.identity <- cell.identity %>% dplyr::select(Sample, Marker, prop.expressed) %>%
      spread(key = Marker, value = prop.expressed) 
  }else if(method == 'mean.expr'){
    cell.identity <- cell.identity %>% dplyr::select(Sample, Marker, mean.expr) %>%
      spread(key = Marker, value = mean.expr) 
  }
  
  cell.identity <- cell.identity %>% mutate(cell.cycle = ifelse(G1 > SM , 'G1',
                                                                ifelse(G1 < SM, 'SM', 'NA')),
                                            life.stage = ifelse(Tachyzoite > Bradyzoite, 'Tachy',
                                                                ifelse(Tachyzoite < Bradyzoite, 'Brady',
                                                                       'NA')))
  cell.identity <- cell.identity %>% dplyr::select(Sample, cell.cycle, life.stage) %>% distinct()
  
  return(cell.identity)
}


getCellIdentity.V2 <- function(expr.norm, method = 'mix'){
  cell.identity <- expr.norm %>% ungroup() %>% dplyr::select(Sample, expr, quantile, Marker) %>% distinct() 
  cell.identity <- cell.identity %>% dplyr::filter(!is.na(Marker))
  cell.identity <- cell.identity %>% group_by(Sample, Marker) %>% 
    summarise(mean.expr = mean(expr),
              sd.expr = sd(expr),
              max.quantile = max(quantile), 
              prop.expressed = sum(quantile >= 2) / n(),
              prop.highly.expressed = sum(quantile >= 3) / n(),
              prop.repressed = sum(quantile < 2) / n())
  cell.identity <- cell.identity %>% dplyr::select(Sample, Marker, mean.expr, prop.highly.expressed) %>%
    pivot_wider(names_from = Marker, values_from = c("mean.expr", "prop.highly.expressed")) 
  
  if(method == 'prop.highly.expressed'){
    cell.identity <- cell.identity %>% mutate(cell.cycle = ifelse(prop.highly.expressed_G1 > prop.highly.expressed_SM , 'G1',
                                                                  ifelse(prop.highly.expressed_G1 < prop.highly.expressed_SM, 'SM', 'NA')),
                                              life.stage = ifelse(prop.highly.expressed_Tachyzoite > prop.highly.expressed_Bradyzoite, 'Tachy',
                                                                  ifelse(prop.highly.expressed_Tachyzoite < prop.highly.expressed_Bradyzoite, 'Brady',
                                                                         'NA')))
    
  }else if(method == 'mean.expr'){
    cell.identity <- cell.identity %>% mutate(cell.cycle = ifelse(mean.expr_G1 > mean.expr_SM , 'G1',
                                                                  ifelse(mean.expr_G1 < mean.expr_SM, 'SM', 'NA')),
                                              life.stage = ifelse(mean.expr_Tachyzoite > mean.expr_Bradyzoite, 'Tachy',
                                                                  ifelse(mean.expr_Tachyzoite < mean.expr_Bradyzoite, 'Brady',
                                                                         'NA')))
    
  }else if(method == 'mix'){
    cell.identity <- cell.identity %>% mutate(cell.cycle = ifelse(prop.highly.expressed_G1 > prop.highly.expressed_SM , 'G1',
                                                                  ifelse(prop.highly.expressed_G1 < prop.highly.expressed_SM, 'SM', 'NA')),
                                              life.stage = ifelse(mean.expr_Tachyzoite > mean.expr_Bradyzoite, 'Tachy',
                                                                  ifelse(mean.expr_Tachyzoite < mean.expr_Bradyzoite, 'Brady',
                                                                         'NA')))
  }
  
  cell.identity <- cell.identity %>% dplyr::select(Sample, cell.cycle, life.stage) %>% distinct()
  
  return(cell.identity)
}


getCellIdentity.V3 <- function(expr.norm, method = 'mix'){
  cell.identity <- expr.norm %>% ungroup() %>% dplyr::select(Sample, expr, quantile, Marker) %>% distinct() 
  cell.identity <- cell.identity %>% dplyr::filter(!is.na(Marker))
  cell.identity <- cell.identity %>% group_by(Sample, Marker) %>% 
    summarise(mean.expr = mean(expr),
              sd.expr = sd(expr),
              max.quantile = max(quantile), 
              prop.expressed = sum(quantile >= 2) / n(),
              prop.highly.expressed = sum(quantile >= 3) / n(),
              prop.repressed = sum(quantile < 2) / n())
  cell.identity <- cell.identity %>% dplyr::select(Sample, Marker, mean.expr, prop.highly.expressed) %>%
    pivot_wider(names_from = Marker, values_from = c("mean.expr", "prop.highly.expressed")) 
  
  if(method == 'prop.highly.expressed'){
    cell.identity <- cell.identity %>% mutate(cell.cycle = ifelse(prop.highly.expressed_G1 > prop.highly.expressed_SM , 'G1',
                                                                  ifelse(prop.highly.expressed_G1 < prop.highly.expressed_SM, 'SM', 'NA')),
                                              life.stage = ifelse(prop.highly.expressed_Tachyzoite > prop.highly.expressed_Bradyzoite, 'Tachy',
                                                                  ifelse(prop.highly.expressed_Tachyzoite < prop.highly.expressed_Bradyzoite, 'Brady',
                                                                         'NA')))
    
  }else if(method == 'mean.expr'){
    cell.identity <- cell.identity %>% mutate(cell.cycle = ifelse(mean.expr_G1 > mean.expr_SM , 'G1',
                                                                  ifelse(mean.expr_G1 < mean.expr_SM, 'SM', 'NA')),
                                              life.stage = ifelse(mean.expr_Tachyzoite > mean.expr_Bradyzoite, 'Tachy',
                                                                  ifelse(mean.expr_Tachyzoite < mean.expr_Bradyzoite, 'Brady',
                                                                         'NA')))
    
  }else if(method == 'mix'){
    cell.identity <- cell.identity %>% mutate(cell.cycle = ifelse(prop.highly.expressed_G1 > prop.highly.expressed_SM , 'G1',
                                                                  ifelse(prop.highly.expressed_G1 < prop.highly.expressed_SM, 'SM', 'NA')),
                                              life.stage = ifelse(mean.expr_Tachyzoite > mean.expr_Bradyzoite, 'Tachy',
                                                                  ifelse(mean.expr_Tachyzoite < mean.expr_Bradyzoite, 'Brady',
                                                                         'NA')))
  }
  
  cell.identity <- cell.identity %>% dplyr::select(Sample, cell.cycle, life.stage) %>% distinct()
  
  return(cell.identity)
}

getUMAP <- function(S.O, cell.identity, title = title, out.file ){
  umap <- S.O[['umap']]@cell.embeddings
  umap <- as.data.frame(umap) %>% mutate(Sample = rownames(umap)) 
  umap <- inner_join(cell.identity, umap, by = 'Sample')
  umap <- umap %>% ungroup() %>% mutate_at(vars(-contains('UMAp')), funs(factor))
  
  p1 <- ggplot(umap, aes(x=UMAP_1,y=UMAP_2)) + 
    geom_point(aes(
      fill = life.stage), shape=21, size = 1) + theme_bw()
  p2 <- ggplot(umap, aes(x=UMAP_1,y=UMAP_2)) + 
    geom_point(aes(
      fill = cell.cycle), shape=21, size = 1) + theme_bw()
  
  p <- grid.arrange(p1, p2, ncol=2, top = textGrob(title))
  ggsave(out.file,
         plot = p, # or give ggplot object name as in myPlot,
         width = 8, height = 5, 
         units = "in", # other options c("in", "cm", "mm"), 
         dpi = 300)
  p
  return(umap)
}

getPCA <- function(S.O, cell.identity, title = title, out.file ){
  pc <- S.O[['pca']]@cell.embeddings
  pc <- as.data.frame(pc) %>% mutate(Sample = rownames(pc)) 
  pc <- inner_join(cell.identity, pc, by = 'Sample')
  pc <- pc %>% ungroup() %>% mutate_at(vars(-contains('PC')), funs(factor))
  
  p1 <- ggplot(pc, aes(x=PC_1,y=PC_2)) + 
    geom_point(aes(
      fill = life.stage), shape=21, size = 1) + theme_bw()
  p2 <- ggplot(pc, aes(x=PC_1,y=PC_2)) + 
    geom_point(aes(
      fill = cell.cycle), shape=21, size = 1) + theme_bw()
  
  p <- grid.arrange(p1, p2, ncol=2, top = textGrob(title))
  ggsave(out.file,
         plot = p, # or give ggplot object name as in myPlot,
         width = 8, height = 5, 
         units = "in", # other options c("in", "cm", "mm"), 
         dpi = 300)
  p
  return(umap)
}

cell.identity.WT.D3 <- getCellIdentity.V2(expr.norm.WT.D3, method = 'mix')
cell.identity.WT.D3 <- getCellIdentity.V2(expr.norm.WT.D3, method = 'mean.expr')
cell.identity.WT.pH <- getCellIdentity(expr.norm.WT.pH, method = 'mean.expr')
cell.identity.KO.D3 <- getCellIdentity(expr.norm.KO.D3, method = 'mean.expr')
cell.identity.KO.pH <- getCellIdentity(expr.norm.KO.pH, method = 'mean.expr')
cell.identity.abundant.markers <- getCellIdentity(expr.norm.abundant.markers, method = 'prop.highly.expressed')
cell.identity.abundant.markers <- getCellIdentity(expr.norm.abundant.markers, method = 'mean.expr')



cell.identity.abundant.markers.WT <- getCellIdentity.V2(expr.norm.abundant.markers.WT, method = 'mix')
cell.identity.abundant.markers.WT.D3 <- getCellIdentity.V2(expr.norm.abundant.markers.WT.D3, method = 'mix')
cell.identity.abundant.markers.WT.pH <- getCellIdentity.V2(expr.norm.abundant.markers.WT.pH, method = 'mix')
cell.identity.abundant.markers <- getCellIdentity.V2(expr.norm.abundant.markers, method = 'mix')
cell.identity.abundant.markers.D3 <- getCellIdentity.V2(expr.norm.abundant.markers.D3, method = 'mix')


umap <- getUMAP(S.O.WT.pH, cell.identity.WT.pH, title = 'WT stressed', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/WT_stressed.pdf')
umap <- getUMAP(S.O.WT.D3, cell.identity.WT.D3, title = 'WT unstressed', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/WT_unstressed.pdf')
umap <- getUMAP(S.O.KO.pH, cell.identity.KO.pH, title = 'KO stressed', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/KO_stressed.pdf')
umap <- getUMAP(S.O.KO.D3, cell.identity.KO.D3, title = 'KO unstressed', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/KO_unstressed.pdf')
umap <- getUMAP(S.O.abundant, cell.identity.abundant, title = 'abundant', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/abundant.pdf')
umap <- getUMAP(S.O.abundant.markers.D3, cell.identity.abundant.markers.D3, title = 'abundant', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/abundant_unstressed.pdf')


umap <- getUMAP(S.O.abundant.markers.WT, cell.identity.abundant.markers.WT, title = 'abundant WT', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/abundant_WT.pdf')
umap <- getUMAP(S.O.abundant.markers.WT.D3, cell.identity.abundant.markers.WT.D3, title = 'abundant WT Unstressed', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/abundant_WT_unstressed.pdf')
umap <- getUMAP(S.O.abundant.markers.WT.pH, cell.identity.abundant.markers.WT.pH, title = 'abundant WT Stressed', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/abundant_WT_stressed.pdf')
umap <- getUMAP(S.O.abundant.markers, cell.identity.abundant.markers, title = 'abundant', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/abundant.pdf')


pc <- getPCA(S.O.abundant.markers.WT.pH, cell.identity.abundant.markers.WT.pH, title = 'abundant WT Stressed', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/abundant_WT_stressed_PCA.pdf')
pc <- getPCA(S.O.abundant.markers.WT.D3, cell.identity.abundant.markers.WT.D3, title = 'abundant WT Unstressed', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/abundant_WT_unstressed_PCA.pdf')
pc <- getPCA(S.O.abundant.markers.WT, cell.identity.abundant.markers.WT, title = 'abundant WT', out.file = '~/work/ToxoPlasmaGondiiR/Output/scFigures/abundant_WT_PCA.pdf')



getMicronemeExpr <- function(expr.norm, cell.identity, microneme.markers){
  expr.norm <- left_join(expr.norm, cell.identity, by = 'Sample')
  microneme.markers.expr <- left_join(microneme.markers, expr.norm, by = c('Gene.ID' = 'GeneID'))
  microneme.markers.expr$cell.cycle[microneme.markers.expr$cell.cycle != 'G1'] = 'NotG1'
  microneme.markers.expr <- microneme.markers.expr %>% 
    dplyr::select(Gene.ID, Name, Sample, expr, cell.cycle)
  microneme.markers.expr <- microneme.markers.expr %>% dplyr::filter(!is.na(expr))
  
  ## Calculate quantiles for each phase
  microneme.markers.expr <- microneme.markers.expr %>% group_by(Gene.ID, cell.cycle) %>% 
    mutate(quantile = ifelse(expr <= quantile(expr)[2], 1, 
                             ifelse(expr <= quantile(expr)[3], 2, 
                                    ifelse(expr <= quantile(expr)[4], 3,4))))
  microneme.markers.expr.summary <- microneme.markers.expr %>%  group_by(Gene.ID) %>% 
    summarise(Name = unique(Name)[1], 
              mean_expr_G1 = mean(expr[cell.cycle == 'G1']),
              sd_expr_G1 = sd(expr[cell.cycle == 'G1']),
              mean_expr_NotG1 = mean(expr[cell.cycle == 'NotG1']),
              sd_expr_NotG1 = sd(expr[cell.cycle == 'NotG1']),
              FC = mean(expr[cell.cycle == 'G1']) / mean(expr[cell.cycle == 'NotG1']),
              pvalue = t.test(expr[cell.cycle == 'G1'], expr[cell.cycle == 'NotG1'])$p.value,
              mean_quantile_G1 = mean(quantile[cell.cycle == 'G1']),
              mean_quantile_NotG1 = mean(quantile[cell.cycle == 'NotG1'])) %>%
    arrange(pvalue)
  
  
  microneme.markers.expr.summary <- microneme.markers.expr.summary %>% na.omit()
  
  #microneme.markers.expr.summary <- microneme.markers.expr.summary  %>% 
  #  pivot_wider(names_from = "cell.cycle", values_from = c("mean_expr", "sd_expr") )
  
  return(microneme.markers.expr.summary)
}


microneme.markers.expr.summary.WT.D3 <- getMicronemeExpr(expr.norm.WT.D3, cell.identity.WT.D3, microneme.markers)
microneme.markers.expr.summary.WT.pH <- getMicronemeExpr(expr.norm.WT.pH, cell.identity.WT.pH, microneme.markers)
microneme.markers.expr.summary.KO.D3 <- getMicronemeExpr(expr.norm.KO.D3, cell.identity.KO.D3, microneme.markers)
microneme.markers.expr.summary.KO.pH <- getMicronemeExpr(expr.norm.KO.pH, cell.identity.KO.pH, microneme.markers)

## Look at samples in G1 phase


write.xlsx(microneme.markers.expr.summary.WT.D3, '~/work/ToxoPlasmaGondiiR/Output/Microneme/microneme_expression_WT_unstressed.xlsx')
write.xlsx(microneme.markers.expr.summary.WT.pH, '~/work/ToxoPlasmaGondiiR/Output/Microneme/microneme_expression_WT_stressed.xlsx')
write.xlsx(microneme.markers.expr.summary.KO.D3, '~/work/ToxoPlasmaGondiiR/Output/Microneme/microneme_expression_KO_unstressed.xlsx')
write.xlsx(microneme.markers.expr.summary.KO.pH, '~/work/ToxoPlasmaGondiiR/Output/Microneme/microneme_expression_KO_stressed.xlsx')


