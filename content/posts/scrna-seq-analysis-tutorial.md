---
title: "单细胞RNA测序数据分析完整流程"
date: 2024-12-05
categories: ["数据分析"]
tags: ["单细胞", "转录组", "R", "Seurat"]
author: "张义博"
description: "从原始数据到细胞类型注释的完整单细胞RNA测序分析教程"
image: uploads/photos/bricks2.jpg
featured: true
---

# 单细胞RNA测序数据分析完整流程

单细胞RNA测序（scRNA-seq）技术已成为研究细胞异质性和发育过程的重要工具。本文将介绍使用Seurat进行scRNA-seq数据分析的完整流程。

## 1. 环境准备

```r
# 安装必要的R包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Seurat")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scater")

# 加载包
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
```

## 2. 数据导入和预处理

### 2.1 读取10X数据

```r
# 读取10X Genomics数据
data_dir <- "path/to/filtered_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir)

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(
  counts = data, 
  project = "scRNA_analysis",
  min.cells = 3,     # 基因至少在3个细胞中表达
  min.features = 200 # 细胞至少表达200个基因
)
```

### 2.2 质量控制指标计算

```r
# 计算线粒体基因比例
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# 计算核糖体基因比例
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")

# 可视化质控指标
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

## 3. 数据过滤

### 3.1 设置过滤条件

```r
# 基于质控指标过滤细胞
seurat_obj <- subset(seurat_obj, subset = 
  nFeature_RNA > 200 & 
  nFeature_RNA < 5000 & 
  percent.mt < 20
)

print(paste("过滤后细胞数:", ncol(seurat_obj)))
print(paste("过滤后基因数:", nrow(seurat_obj)))
```

### 3.2 过滤高变基因

```r
# 识别高变基因
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 可视化高变基因
top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

## 4. 数据标准化和降维

### 4.1 数据标准化

```r
# 标准化数据
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

# 或者只对高变基因进行标准化（更快）
# seurat_obj <- ScaleData(seurat_obj)
```

### 4.2 主成分分析

```r
# 进行PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# 可视化PCA结果
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(seurat_obj, reduction = "pca")
```

### 4.3 确定主成分数量

```r
# 绘制ElbowPlot确定主成分数量
ElbowPlot(seurat_obj, ndims = 50)

# 基于JackStraw方法（可选，计算时间较长）
# seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
# seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
# JackStrawPlot(seurat_obj, dims = 1:15)
```

## 5. 细胞聚类

### 5.1 构建邻接图

```r
# 构建KNN图
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# 进行聚类
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# 查看聚类结果
head(Idents(seurat_obj), 5)
```

### 5.2 UMAP降维可视化

```r
# 运行UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# 可视化聚类结果
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

## 6. 寻找标志基因

### 6.1 寻找每个cluster的标志基因

```r
# 寻找所有cluster的标志基因
cluster.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 显示top基因
cluster.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
```

### 6.2 可视化标志基因

```r
# 热图展示top基因
top5 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(seurat_obj, features = top5$gene) + NoLegend()

# 特征图展示特定基因
FeaturePlot(seurat_obj, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A"))
```

## 7. 细胞类型注释

### 7.1 基于已知标志基因注释

```r
# 定义细胞类型标志基因
cell_type_markers <- list(
  "T cells" = c("CD3E", "CD3D"),
  "B cells" = c("MS4A1", "CD79A"),
  "NK cells" = c("GNLY", "NKG7"),
  "Monocytes" = c("CD14", "LYZ"),
  "Dendritic cells" = c("FCER1A", "CST3")
)

# 手动注释细胞类型
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

# 可视化注释结果
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

### 7.2 使用SingleR自动注释

```r
# 使用SingleR进行自动注释
library(SingleR)
library(celldex)

# 获取参考数据集
ref <- HumanPrimaryCellAtlasData()

# 提取表达矩阵
sce <- as.SingleCellExperiment(seurat_obj)

# 进行注释
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)

# 添加注释结果到Seurat对象
seurat_obj$SingleR_labels <- pred$labels

# 可视化自动注释结果
DimPlot(seurat_obj, reduction = "umap", group.by = "SingleR_labels", label = TRUE)
```

## 8. 差异表达分析

### 8.1 不同细胞类型间的差异分析

```r
# 比较两个细胞类型
tcell_vs_bcell <- FindMarkers(seurat_obj, ident.1 = "T cells", ident.2 = "B cells")
head(tcell_vs_bcell, n = 10)

# 火山图可视化
library(EnhancedVolcano)
EnhancedVolcano(tcell_vs_bcell,
  lab = rownames(tcell_vs_bcell),
  x = 'avg_log2FC',
  y = 'p_val_adj')
```

### 8.2 条件间的差异分析

```r
# 如果有多个条件（如对照组vs处理组）
Idents(seurat_obj) <- "condition"  # 假设有condition这个metadata

# 在特定细胞类型中比较条件
tcell_subset <- subset(seurat_obj, idents = "T cells")
condition_markers <- FindMarkers(tcell_subset, ident.1 = "treatment", ident.2 = "control")
```

## 9. 轨迹分析

### 9.1 使用Monocle3进行轨迹分析

```r
library(monocle3)

# 转换Seurat对象为CDS对象
cds <- as.cell_data_set(seurat_obj)

# 预处理
cds <- preprocess_cds(cds, num_dim = 100)

# 降维
cds <- reduce_dimension(cds)

# 聚类
cds <- cluster_cells(cds)

# 学习轨迹
cds <- learn_graph(cds)

# 排序细胞
cds <- order_cells(cds)

# 可视化轨迹
plot_cells(cds, 
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
```

## 10. 结果输出

### 10.1 保存分析结果

```r
# 保存Seurat对象
saveRDS(seurat_obj, file = "analyzed_seurat_object.rds")

# 导出细胞注释结果
write.csv(seurat_obj@meta.data, "cell_annotations.csv")

# 导出标志基因
write.csv(cluster.markers, "cluster_markers.csv")
```

### 10.2 生成报告图表

```r
# 创建综合图表
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type")
p3 <- FeaturePlot(seurat_obj, features = "nFeature_RNA")
p4 <- VlnPlot(seurat_obj, features = "percent.mt", group.by = "cell_type")

# 组合图表
combined_plot <- (p1 + p2) / (p3 + p4)
ggsave("scRNA_analysis_summary.pdf", combined_plot, width = 12, height = 10)
```

## 总结

本文介绍了使用Seurat进行单细胞RNA测序数据分析的完整流程，包括：

1. **数据预处理**：质量控制和过滤
2. **标准化降维**：数据标准化和主成分分析
3. **聚类分析**：细胞聚类和可视化
4. **细胞注释**：标志基因识别和细胞类型注释
5. **下游分析**：差异表达和轨迹分析

这个流程可以根据具体的研究问题和数据特点进行调整和优化。

## 参考资料

1. [Seurat官方教程](https://satijalab.org/seurat/)
2. [单细胞数据分析最佳实践](https://www.sc-best-practices.org/)
3. [Orchestrating Single-Cell Analysis with Bioconductor](https://osca.bioconductor.org/)
