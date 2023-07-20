setwd("/Users/liuzhe/Desktop/cityu/COVID/02_ImmuneInfiltraction")
remove(list=ls())

library("xCell")
#https://www.yunbios.net/cn/xCell.html
#该方法适用于基因表达谱和传统RNA-seq数据，但不包括单细胞数据（推荐singleR，同样由该团队开发）。

#输入数据应该是RPKM/FPKM/TPM/RSEM 而不是 row read counts
#https://xcell.ucsf.edu/
#To use xCell simply upload human gene expression data file in tab delimited text format or csv (up to 1Gb). The expression matrix should be a matrix with genes in rows and samples in columns. The rownames should be gene symbols. If the data contains non-unique gene symbols, rows with same gene symbols will be averaged. xCell uses the expression levels ranking and not the actual values, thus normalization does not have an effect, however normalizing to gene length (RPKM/FPKM/TPM/RSEM) is required.

exp<-read.csv("/Users/liuzhe/Desktop/cityu/COVID/01_DEGs/discoverydataset_GSE157103_GSE152641.csv", row.names = 1)
# GSE157103和GSE152641都是bulk RNA-seq数据
dim(exp)
#[1] 19184   212
exp[1:5,1:5]
#62列是COVID19 24列是healthy controls 100列是COVID9 26列是healthy controls

#counts转化为FPKM
library("biomaRt")
ensembl <- useMart("ensembl") 
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
test <- getBM(attributes=c('ensembl_gene_id', 'start_position',
                           'end_position','ensembl_transcript_id',
                           'transcript_length','hgnc_symbol'),mart = ensembl)
test <- test[order(test$transcript_length,decreasing = T),]
g <- test[!duplicated(test$hgnc_symbol),]
g <- g[,c(5,6)]
head(g)
dim(g)
summary(g)

exp$hgnc_symbol<-rownames(exp)
counts_geneId <- merge(exp,g,by = "hgnc_symbol")
rownames(counts_geneId)<-counts_geneId$hgnc_symbol
counts_geneId$hgnc_symbol<-NULL
kb<-counts_geneId$transcript_length/1000
countdata<-counts_geneId[,1:212]
rpk<-countdata/kb
fpkm<-t(t(rpk)/colSums(countdata)*10^6)
dim(fpkm)
#[1] 19060   212
fpkm[1:5,1:4]
write.csv(fpkm, file = "FPKM_discoverydataset_GSE157103_GSE152641.csv", quote = F)
fpkm<-as.data.frame(fpkm)

## 如果是RNA-seq数据，则
xCell_RNAseq<- xCellAnalysis(fpkm,rnaseq = T)
dim(xCell_RNAseq)
#[1]  67 212
xCell_RNAseq[1:5,1:5]
xCell_RNAseq<-as.data.frame(xCell_RNAseq)
xCell_RNAseq_resort<-subset(xCell_RNAseq,select = c(1:62, 87:186, 63:86, 187:212))
dim(xCell_RNAseq_resort)
#[1]  67 212
write.csv(xCell_RNAseq_resort, "Cell_score.csv", quote = F)
#前162个样本属于COVID19组；后50个样本属于healthy controls组

#https://blog.csdn.net/u011794835/article/details/124943501
#https://www.jianshu.com/p/6f558be1b03a
library(ggpubr)
library(rstatix)
#xCell_RNAseq_resort<-read.csv("Cell_score.csv", header = T, row.names = 1)
inputdata<-t(xCell_RNAseq_resort)
inputdata<-as.data.frame(inputdata)
inputdata$group<-c(rep("COVID19",162),rep("HealthControl",50))
inputdata$group<-as.factor(inputdata$group)
inputdata$sampleName<-rownames(inputdata)

#数据变形 melt()
library(reshape2)
library(knitr)
mydata<-melt(inputdata, id.vars = c("sampleName","group"), variable.name = "cellType", value.name = "xCellScore")

myt_test<-t_test(group_by(mydata, cellType), xCellScore~group)
myt_test<-add_significance(myt_test, "p")
my_t.test<-add_xy_position(myt_test, x = "cellType", dodge =  0.8)
results_sig<-as.matrix(my_t.test)
write.csv(results_sig, "difference_sig_levels_t_test.csv", row.names = F, quote = F)
celltypes_uniq<-unique(my_t.test$cellType)
#write.csv(celltypes_uniq, "celltypes_uniq.csv", row.names = F, quote = F)

mywilcox_test<-wilcox_test(group_by(mydata, cellType), xCellScore~group)
mywilcox_test<-add_significance(mywilcox_test, "p")
my_mywilcox_test<-add_xy_position(mywilcox_test, x = "cellType", dodge =  0.8)
results_sig<-as.matrix(my_mywilcox_test)
write.csv(results_sig, "difference_sig_levels_wilcox_test.csv", row.names = F, quote = F)


#一共64种细胞类型，其中包括了14个基质细胞； 22个免疫细胞
#paper Figure1
#Lymphoids
#Stem cells
#Myeloids
#Stromal cells
#Others

lymphoids<-c("B-cells", "CD4+ memory T-cells", "CD4+ naive T-cells", "CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem",
             "CD8+ naive T-cells", "CD8+ T-cells", "CD8+ Tcm", "CD8+ Tem", "Class-switched memory B-cells",
             "Memory B-cells", "naive B-cells", "NK cells", "NKT", "Plasma cells", "pro B-cells", "Tgd cells",
             "Th1 cells", "Th2 cells", "Tregs")
length(lymphoids)
#[1] 21
stemcells<-c("CLP", "CMP", "Erythrocytes", "GMP", "HSC", "Megakaryocytes", "MEP", "MPP", "Platelets")
length(stemcells)
#[1] 9
myeloids<-c("aDC", "Basophils", "cDC", "DC", "Eosinophils", "iDC", "Macrophages", "Macrophages M1", "Macrophages M2",
            "Mast cells", "Monocytes", "Neutrophils", "pDC")
length(myeloids)
#[1] 13
stromalcells<-c("Adipocytes", "Chondrocytes", "Endothelial cells", "Fibroblasts", "ly Endothelial cells", "MSC",
                "mv Endothelial cells", "Myocytes", "Osteoblast", "Pericytes", "Preadipocytes", "Skeletal muscle",
                "Smooth muscle")
length(stromalcells)
#[1] 13
others<-c("Astrocytes", "Epithelial cells", "Hepatocytes", "Keratinocytes", "Melanocytes", "Mesangial cells", "Neurons",
          "Sebocytes")
length(others)
#[1] 8

library(tidyverse)
##############Class1: Lymphoids###########################################
mydata_selected <- mydata %>% filter(cellType %in% lymphoids)
#myt_test<-t_test(group_by(mydata_selected, cellType), xCellScore~group)
#myt_test<-add_significance(myt_test, "p")
#my_t.test<-add_xy_position(myt_test, x = "cellType", dodge =  0.8)
pdf("EnrichmentScore_lymphoids_sig.pdf", height = 6, width = 8)
p <- ggboxplot(mydata_selected, x = "cellType", y = "xCellScore",fill = "group", color = 'group',
               palette = c("#00AFBB", "#E7B800"), add = "boxplot", add.params = list(fill = "white")) +
  labs(x = "Cell type", y = "Cell score from xCell", fill = "group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = group), label = "p.signif", paired = FALSE)
dev.off()
pdf("EnrichmentScore_lymphoids_pval.pdf")
p <- ggboxplot(mydata_selected, x = "cellType", y = "xCellScore",fill = "group", color = 'group',
               palette = c("#00AFBB", "#E7B800"), add = "boxplot", add.params = list(fill = "white")) +
  labs(x = "Cell type", y = "Cell score from xCell", fill = "group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = group), label = "p.format", paired = FALSE)
dev.off()

##############Class2: Stem Cells##########################################
mydata_selected <- mydata %>% filter(cellType %in% stemcells)
pdf("EnrichmentScore_stemcells_sig.pdf", height = 5, width = 4)
p <- ggboxplot(mydata_selected, x = "cellType", y = "xCellScore",fill = "group", color = 'group',
               palette = c("#00AFBB", "#E7B800"), add = "boxplot", add.params = list(fill = "white")) +
  labs(x = "Cell type", y = "Cell score from xCell", fill = "group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = group), label = "p.signif", paired = FALSE)
dev.off()
pdf("EnrichmentScore_stemcells_pval.pdf")
p <- ggboxplot(mydata_selected, x = "cellType", y = "xCellScore",fill = "group", color = 'group',
               palette = c("#00AFBB", "#E7B800"), add = "boxplot", add.params = list(fill = "white")) +
  labs(x = "Cell type", y = "Cell score from xCell", fill = "group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = group), label = "p.format", paired = FALSE)
dev.off()

##############Class3: Myeloids############################################
mydata_selected <- mydata %>% filter(cellType %in% myeloids)
pdf("EnrichmentScore_myeloids_sig.pdf", height = 5, width = 4)
p <- ggboxplot(mydata_selected, x = "cellType", y = "xCellScore",fill = "group", color = 'group',
               palette = c("#00AFBB", "#E7B800"), add = "boxplot", add.params = list(fill = "white")) +
  labs(x = "Cell type", y = "Cell score from xCell", fill = "group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = group), label = "p.signif", paired = FALSE)
dev.off()
pdf("EnrichmentScore_myeloids_pval.pdf")
p <- ggboxplot(mydata_selected, x = "cellType", y = "xCellScore",fill = "group", color = 'group',
               palette = c("#00AFBB", "#E7B800"), add = "boxplot", add.params = list(fill = "white")) +
  labs(x = "Cell type", y = "Cell score from xCell", fill = "group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = group), label = "p.format", paired = FALSE)
dev.off()

##############Class4: Stromal Cells#######################################
mydata_selected <- mydata %>% filter(cellType %in% stromalcells)
pdf("EnrichmentScore_stromalcells_sig.pdf", height = 5, width = 4)
p <- ggboxplot(mydata_selected, x = "cellType", y = "xCellScore",fill = "group", color = 'group',
               palette = c("#00AFBB", "#E7B800"), add = "boxplot", add.params = list(fill = "white")) +
  labs(x = "Cell type", y = "Cell score from xCell", fill = "group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = group), label = "p.signif", paired = FALSE)
dev.off()
pdf("EnrichmentScore_stromalcells_pval.pdf")
p <- ggboxplot(mydata_selected, x = "cellType", y = "xCellScore",fill = "group", color = 'group',
               palette = c("#00AFBB", "#E7B800"), add = "boxplot", add.params = list(fill = "white")) +
  labs(x = "Cell type", y = "Cell score from xCell", fill = "group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = group), label = "p.format", paired = FALSE)
dev.off()

##############Class5: Others##############################################
mydata_selected <- mydata %>% filter(cellType %in% others)
pdf("EnrichmentScore_others_sig.pdf", height = 5, width = 4)
p <- ggboxplot(mydata_selected, x = "cellType", y = "xCellScore",fill = "group", color = 'group',
               palette = c("#00AFBB", "#E7B800"), add = "boxplot", add.params = list(fill = "white")) +
  labs(x = "Cell type", y = "Cell score from xCell", fill = "group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = group), label = "p.signif", paired = FALSE)
dev.off()
pdf("EnrichmentScore_others_pval.pdf")
p <- ggboxplot(mydata_selected, x = "cellType", y = "xCellScore",fill = "group", color = 'group',
               palette = c("#00AFBB", "#E7B800"), add = "boxplot", add.params = list(fill = "white")) +
  labs(x = "Cell type", y = "Cell score from xCell", fill = "group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = group), label = "p.format", method = "t.test", paired = FALSE)
dev.off()

##统计学方法：wilcox.test
#ns: p > 0.05
#*: p <= 0.05
#**: p <= 0.01
#***: p <= 0.001
#****: p <= 0.0001
#cellType        #no_of_celltypes    #not_significant
# Lymphoids                21              4
# Stem Cells               9               5
# Myeloids                 13              6
# Stromal Cells            13              8
# Others                    8              6




