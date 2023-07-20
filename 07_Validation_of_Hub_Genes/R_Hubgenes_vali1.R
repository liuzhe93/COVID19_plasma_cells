setwd("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/07_Hubgenes")
remove(list=ls())

hubgenes<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/06_PPI/top20_detail.csv")
hub<-hubgenes$X


#################################validation dataset1##################################
counts_GSE152418<-read.table("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/DataCollection/COVID19_bulkRNA/GSE152418/GSE152418_p20047_Study1_RawCounts.txt",
                             row.names = 1, header = T, sep = "\t")
dim(counts_GSE152418)
#[1] 60683    34
clinical_GSE152418<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/DataCollection/COVID19_bulkRNA/GSE152418/clinical.csv",
                             header = T)

library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DESeq2)

counts_GSE152418$geneid <- rownames(counts_GSE152418)
# transform id  

map_dt <- bitr(counts_GSE152418$geneid, fromType = "ENSEMBL",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
dt_merge <- merge(map_dt,counts_GSE152418, by.y = "geneid", by.x = "ENSEMBL")
dt_merge_uniq<-aggregate(dt_merge[,3:36], by = list(type = dt_merge$SYMBOL), mean)
rownames(dt_merge_uniq) <- dt_merge_uniq$type
dt_merge_uniq$type <- NULL
dim(dt_merge_uniq)
#[1] 35432    34
write.csv(dt_merge_uniq, "validationdataset_GSE152418.csv", quote = F)


condition <- clinical_GSE152418$severity
sampleTable <- data.frame(condition=condition)
#样本表达数据框列名需要与sampleTable一致。
row.names(sampleTable) <- colnames(dt_merge_uniq) 
data = apply(dt_merge_uniq, 2, as.integer) ## DESeq2分析需要是整数
row.names(data) <- row.names(dt_merge_uniq)

dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleTable, design = ~condition)
dds$condition<- relevel(dds$condition, ref = "Healthy") 
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- vst(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
valida_1<-normal_gset[hub,]
val_dat1<-t(valida_1)
val_dat1<-as.data.frame(val_dat1)
sel_dat1<-cbind(val_dat1, clinical_GSE152418$severity)
dim(sel_dat1)
#[1] 34 21
colnames(sel_dat1)[21]<-"group"

library("ggplot2")
library("ggpubr")
library("ggstatsplot")
library("ggsci")

pdf("vali_1_CCNA2_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("CCNA2", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = CCNA2, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_TOP2A_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("TOP2A", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = TOP2A, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()


pdf("vali_1_CCNB2_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("CCNB2", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = CCNB2, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_TPX2_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("TPX2", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = TPX2, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()


pdf("vali_1_BUB1_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("BUB1", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = BUB1, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_CDK1_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("CDK1", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = CDK1, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_KIF2C_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("KIF2C", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = KIF2C, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_RRM2_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("RRM2", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = RRM2, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_AURKB_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("AURKB", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = AURKB, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()


pdf("vali_1_CDCA8_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("CDCA8", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = CDCA8, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_KIF11_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("KIF11", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = KIF11, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_CCNB1_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("CCNB1", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = CCNB1, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_CDC20_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("CDC20", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = CDC20, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()


pdf("vali_1_UBE2C_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("UBE2C", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = UBE2C, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_KIF20A_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("KIF20A", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = KIF20A, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_PTTG1_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("PTTG1", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = PTTG1, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_TTK_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("TTK", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = TTK, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_DLGAP5_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("DLGAP5", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = DLGAP5, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_BUB1B_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("BUB1B", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = BUB1B, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

pdf("vali_1_PBK_stage.pdf", width = 4, height = 3)
mydata<-subset(sel_dat1, select = c("PBK", "group"))
ggplot(data = mydata, mapping = aes(x = group, y = PBK, fill= group)) +
  geom_violin() + geom_boxplot(width=0.2)+
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()



