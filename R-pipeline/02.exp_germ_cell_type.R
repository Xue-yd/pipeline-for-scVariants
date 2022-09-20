## calculate the mean expression level of every gene in every cell type

library(Seurat)

fca_germ <- readRDS("fca_germ_seurat.rds")

fca_germ@active.ident <- as.factor(fca_germ$my_annotation)

## 1. mean UMI count
fca_germ_seurat <- NormalizeData(fca_germ_seurat,scale.factor = 10000)
UMI_mean_exp <- AverageExpression(fca_germ_seurat,
                                  group.by = "my_annotation",
                                  slot = "data")
write.table(UMI_mean_exp$RNA,file = "mean_UMI_count.txt",
            sep = "\t")

# 2. scaled expression
all.genes <- rownames(fca_germ_seurat)
fca_germ_seurat <- ScaleData(fca_germ_seurat, features = all.genes)

scaled_exp <- AverageExpression(fca_germ_seurat,
                               group.by = "my_annotation",
                               slot = "scale.data")
write.table(scaled_exp$RNA,file = "mean_scaled_expression.txt",
            sep = "\t")
