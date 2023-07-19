setwd("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/04_scRNA/")
remove(list=ls())

## construct RNA Seurat object
dir.10x = 'download/'
genes = read.delim(paste0(dir.10x, 'lung_geneNames_upload.csv'), header = F)
barcodes = read.delim(paste0(dir.10x, 'lung_cellNames.csv'), header = F)
mtx = Matrix::readMM(paste0(dir.10x, 'gene_sorted-lung_expression_data.mtx'))
mtx = as.matrix(mtx)
mtx = as(mtx, 'dgCMatrix') #将稀疏矩阵转换为Seurat默认的dgC格式
colnames(mtx) = barcodes$V1
rownames(mtx) = genes$V1

library("Seurat")
packageVersion("Seurat")
#[1] ‘4.2.0’
covid19 = CreateSeuratObject(counts = mtx, assay = 'RNA', meta.data = barcodes)
covid19
#An object of class Seurat 
#34546 features across 116314 samples within 1 assay 
#Active assay: RNA (34546 features, 0 variable features)
save(covid19, file = "covid19.RData")

load("covid19.RData")
metadata_lung<-read.csv("download/lung_metaData.txt", header = T, sep = "\t")
dim(metadata_lung)
#[1] 116314     23
metadata_lung_sel<-metadata_lung[2:116314,]
overlap<-intersect(rownames(covid19@meta.data),metadata_lung_sel$NAME)
length(overlap)
#[1] 116313

covid19[["CellName"]] <- colnames(covid19)
covid19_subset<-subset(covid19, subset = CellName %in% overlap)


covid19_subset <- AddMetaData(object = covid19_subset, metadata = metadata_lung_sel$biosample_id, col.name = "biosample_id") 
covid19_subset <- AddMetaData(object = covid19_subset, metadata = metadata_lung_sel$disease__ontology_label, col.name = "disease__ontology_label") 
covid19_subset <- AddMetaData(object = covid19_subset, metadata = metadata_lung_sel$group, col.name = "group") 
covid19_subset <- AddMetaData(object = covid19_subset, metadata = metadata_lung_sel$cell_type_main, col.name = "cell_type_main") 
covid19_subset <- AddMetaData(object = covid19_subset, metadata = metadata_lung_sel$cell_type_intermediate, col.name = "cell_type_intermediate") 
covid19_subset <- AddMetaData(object = covid19_subset, metadata = metadata_lung_sel$cell_type_fine, col.name = "cell_type_fine") 
covid19_subset <- AddMetaData(object = covid19_subset, metadata = metadata_lung_sel$age, col.name = "age") 
covid19_subset <- AddMetaData(object = covid19_subset, metadata = metadata_lung_sel$sex, col.name = "sex")
covid19_subset <- AddMetaData(object = covid19_subset, metadata = metadata_lung_sel$recorded_race, col.name = "recorded_race")
covid19_subset <- AddMetaData(object = covid19_subset, metadata = metadata_lung_sel$intubation_days, col.name = "intubation_days") 
save(covid19_subset, file = "covid19_subset.RData")

Idents(covid19_subset)<-covid19_subset$cell_type_intermediate
save(covid19_subset,file="covid19_subset.anno.RData")


# find all markers of cluster "Plasma cells"
markers <- FindMarkers(covid19_subset, ident.1 = "Plasma cells", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                       test.use = "wilcox")
head(markers, n = 5)
write.table(markers,"covid19_subset_plasma_cells.markers.csv",sep=",",quote=F)
save(markers,file="covid19_subset.markers.RData")


umap_loc<-read.table("download/lung_clusterfile.txt", sep = "\t", header = T)
umap_loc<-umap_loc[-1,]
dim(umap_loc)
#[1] 116313      3


covid19_subset <- AddMetaData(object = covid19_subset, metadata = row.names(covid19_subset@meta.data), col.name = "NAME") 
merged_data<-merge(covid19_subset@meta.data, umap_loc, by = "NAME")
merged_data$X<-as.numeric(merged_data$X)
merged_data$Y<-as.numeric(merged_data$Y)
library("ggplot2")
pdf("UMAP.pdf", width = 11, height = 5)
ggplot(merged_data, aes(x = X, y = Y, colour = cell_type_intermediate)) + geom_point(size = 0.0000000001) + theme_bw() +
#  scale_color_gradientn(values = seq(0,1,0.2),colours = c("blue","cyan","green","yellow","orange","red"))+
  facet_wrap(~group)
dev.off()

prop_celltype<-prop.table(table(Idents(covid19_subset), covid19_subset$group))
write.csv(prop_celltype, 'prop.csv',quote = T,row.names = T)

counts_celltype <- table(Idents(covid19_subset),covid19_subset$group)
write.csv(counts_celltype, 'counts.csv',quote = T,row.names = T)


library("dplyr")
library("reshape2")
library("plyr")

CellRatio<-prop.table(table(Idents(covid19_subset),covid19_subset$biosample_id),margin = 2)
CellRatio
CellRatio<-as.data.frame(CellRatio)
colourCount <- length(unique(CellRatio$Var1))
pdf("stacked_barplot.pdf")
ggplot(CellRatio) + geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = "identity", colour = "#222222") + 
  theme_classic() + labs(x = "BioSample_ID", y = "Cell Ratio") + 
  coord_flip() + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))
dev.off()



cellper <- dcast(CellRatio,Var2~Var1, value.var = "Freq")
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
sample <- unique(covid19_subset@meta.data$biosample_id)
group <- c(rep("Control",7), rep("COVID-19",20))
samples <- data.frame(sample, group)#创建数据框

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group']#R添加列

pplist = list()
sce_groups = unique(covid19_subset@meta.data$cell_type_intermediate)
library("ggpubr")
library("cowplot")
for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))#选择一组数据
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###组间t检验分析
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("Control", "COVID-19") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf("statistics_1_8.pdf")
plot_grid(pplist[['Airway epithelial cells']],pplist[['Macrophages']],pplist[['AT1']],
          pplist[['Smooth muscle']],pplist[['AT2']],pplist[['Fibroblasts']],
          pplist[['Dendritic cells']],pplist[['Cycling NK/T cells']],pplist[['Endothelial cells']])
dev.off()
pdf("statistics_9_18.pdf")
plot_grid(pplist[['Other epithelial cells']], pplist[['Mast cells']],pplist[['Neuronal cells']],
          pplist[['Plasma cells']],pplist[['Monocytes']],pplist[['NK cells']],
          pplist[['CD4+ T cells']],pplist[['B cells']],pplist[['Tregs']])
dev.off()
pdf("statistics_19.pdf")
plot_grid(pplist[['CD8+ T cells']])
dev.off()
