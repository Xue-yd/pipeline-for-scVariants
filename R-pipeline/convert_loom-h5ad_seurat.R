## convert loom/h5ad file to Seurat Object

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(dplyr)
library(tidyverse)
library(purrr)

# 1.testis seurat: all cell type,and reductions information
# 1-1 loom file
ds <- Connect("r_fca_biohub_testis_10x.loom",mode = "r")

fca_testis_mat <- ds[["/matrix"]][,]
fca_testis_mat <- Matrix::Matrix(fca_testis_mat, sparse=T)
fca_testis_mat <- Matrix::t(fca_testis_mat)

fca_testis_cellid <- ds[["/col_attrs/CellID"]][]
fca_testis_geneid <- ds[["/row_attrs/Gene"]][]
colnames(fca_testis_mat) <- fca_testis_cellid
rownames(fca_testis_mat) <- fca_testis_geneid

attrs <- c('CellID', 'ClusterID', 'n_counts', 'n_genes', 'percent_mito','sex',
           'annotation', 'annotation_broad', 'S_annotation_broad_extrapolated')

attrs_df <- map_dfc(attrs, ~ ds[[paste0("col_attrs/", .)]][]) %>% as.data.frame()
colnames(attrs_df) <- attrs
rownames(attrs_df) <- fca_testis_cellid

write.table(attrs_df,file = "attrs.txt",sep = ",",
            row.names = FALSE,quote = F)


# 1-2 h5ad file:add reduction information
Convert("r_fca_biohub_testis_10x.h5ad", dest = "h5Seurat",overwrite = F)
fca_testis_h5 <- LoadH5Seurat("r_fca_biohub_testis_10x.h5seurat")

fca_testis_seurat <- CreateSeuratObject(counts = fca_testis_mat,
                                      meta.data = attrs_df)

fca_testis_seurat@reductions <- fca_testis_h5@reductions
rownames(fca_testis_seurat@reductions$pca@cell.embeddings) <- fca_testis_cellid
rownames(fca_testis_seurat@reductions$tsne@cell.embeddings) <- fca_testis_cellid
rownames(fca_testis_seurat@reductions$umap@cell.embeddings) <- fca_testis_cellid

DimPlot(fca_testis_seurat, 
        reduction = 'tsne', 
        group.by = 'annotation_broad', 
        label = T, 
        label.size = 4) + NoLegend()

write.table(fca_testis_seurat$annotation,file = "cell_ann.txt")
saveRDS(fca_testis_seurat,file = "fca_testis.seurat.rds")

## 2. germ seurat:sub-object

# 2-1.testis matrix
ds <- Connect("r_fca_biohub_testis_10x.loom",mode = "r")
fca_testis_mat <- ds[["/matrix"]][,]

# 2-2.raw counts rename
fca_testis_cellid <- ds[["/col_attrs/CellID"]][]
fca_testis_geneid <- ds[["/row_attrs/Gene"]][]
colnames(fca_testis_mat) <- fca_testis_geneid
rownames(fca_testis_mat) <- fca_testis_cellid
rm(fca_testis_cellid)

# 2-3.extract germ cell raw counts 
germcell <- read.table(file = "attrs_germline.txt",header = T,
                       sep = "\t")
germcell <- germcell$CellID
germcell <- as.array(germcell)
fca_germ_mat <- fca_testis_mat[rownames(fca_testis_mat) %in% germcell,]
rm(fca_testis_mat)

## 2-4.sparse matrix
fca_germ_mat <- Matrix::Matrix(fca_germ_mat, sparse=T)
fca_germ_mat <- Matrix::t(fca_germ_mat)

# 2-5.meta infor 
attrs <- c('CellID', 'ClusterID', 'n_counts', 'n_genes', 'percent_mito','sex',
           'annotation', 'annotation_broad', 'S_annotation_broad_extrapolated')
## 2-5-1
attrs_df <- map_dfc(attrs, ~ ds[[paste0("col_attrs/", .)]][]) %>% as.data.frame()
colnames(attrs_df) <- attrs
rownames(attrs_df) <- fca_testis_cellid

attrs_df <- attrs_df[rownames(attrs_df) %in% germcell,]

## 2-5-2
attrs <- c('CellID', 'ClusterID', 'n_counts', 'n_genes', 'percent_mito','sex',
           'annotation', 'my_annotation','annotation_broad', 'S_annotation_broad_extrapolated')

attrs_df <- read.table(file = "attrs_germline.txt.csv",sep = ",",
                       header = T)
colnames(attrs_df) <- attrs
rownames(attrs_df) <- attrs_df$CellID

fca_germ_seurat <- CreateSeuratObject(counts = fca_germ_mat,
                                        meta.data = attrs_df)

fca_germ_seurat@active.ident <- as.factor(fca_germ_seurat$my_annotation)
saveRDS(fca_germ_seurat,file = "fca_germ_seurat.rds")

