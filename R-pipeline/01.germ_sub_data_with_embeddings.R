library(Seurat)

testisdata <- readRDS("fca_testis.seurat.rds")
mydata <- readRDS("fca_germ_seurat.rds")

germ_data <- subset(x = testisdata, cells = colnames(mydata))
germ_data@active.ident <- mydata@active.ident
germ_data@meta.data$my_annotation <- mydata@active.ident

saveRDS(germ_data,file = "fca_germ_seurat_with_embeddings.rds")
