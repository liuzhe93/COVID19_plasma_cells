setwd("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/08_ROC")
remove(list=ls())

hubgenes<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/06_PPI/top20_detail.csv")
hub<-hubgenes$X

## predict immune response
#ROC: receiver operating characteristic,ROC曲线
#AUC: area under the ROC curve,曲线下面积
#pAUC: partial area under the ROC curve 部分曲线下面积
#CI: confidence interval 可信区间
#SP: specificity 特异度
#SE: sensitivity 灵敏度

require(pROC)




#######################valida2
counts_GSE179627<-read.table("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/DataCollection/COVID19_bulkRNA/GSE179627/GSE179627_gene_reads_count.txt",
                             header = T, row.names = 1, sep = "\t", encoding = 'utf-8', fill = T)
dim(counts_GSE179627)
#[1] 55765    73
clinical_GSE179627<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/DataCollection/COVID19_bulkRNA/GSE179627/clinical.csv",
                             header = T, row.names = 1)
samples<-rownames(clinical_GSE179627)
samples<-gsub("-", "_", samples)
rownames(clinical_GSE179627)<-samples
conditions<-clinical_GSE179627$immune.response
conditions<-gsub(" ", "_", conditions)
conditions<-gsub("-", "_", conditions)
clinical_GSE179627$immune.response<-conditions
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DESeq2)

counts_GSE179627$description<-NULL
counts_GSE179627$gene_biotype<-NULL
counts_GSE179627$X4<-NULL
#样本表达数据框列名需要与sampleTable一致。

srt_sample<-colnames(counts_GSE179627)
clinical_corr<-clinical_GSE179627[srt_sample,]
condition <- clinical_corr$immune.response
sampleTable <- data.frame(condition=condition)
#样本表达数据框列名需要与sampleTable一致。
row.names(sampleTable) <- colnames(counts_GSE179627) 

data = apply(counts_GSE179627, 2, as.integer) ## DESeq2分析需要是整数
row.names(data) <- row.names(counts_GSE179627)

dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleTable, design = ~condition)
dds$condition<- relevel(dds$condition, ref = "uninfected") 
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- vst(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
valida_2<-normal_gset[hub,]

val_dat2<-t(valida_2)
val_dat2<-as.data.frame(val_dat2)
sel_dat2<-cbind(val_dat2, clinical_corr$immune.response)
dim(sel_dat2)
#[1] 70 21
colnames(sel_dat2)[21]<-"group"

library("ggplot2")
library("ggpubr")
library("ggstatsplot")
library("ggsci")

head(sel_dat2)
roc.list <- roc(group ~ CCNA2+TOP2A+CCNB2+TPX2+BUB1+CDK1+KIF2C+RRM2+AURKB+CDCA8+KIF11+CCNB1+CDC20+UBE2C+KIF20A+PTTG1+TTK+DLGAP5+BUB1B+PBK, data = sel_dat2)

####################CCNA2###############################
roc.list$CCNA2
#Call:
#  roc.formula(formula = group ~ CCNA2, data = sel_dat2)

#Data: CCNA2 in 9 controls (group asymptomatic) > 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.5602
ci(roc.list$CCNA2)
#95% CI: 0.2917-0.8286 (DeLong)

####################TOP2A###############################
roc.list$TOP2A
#Call:
#  roc.formula(formula = group ~ TOP2A, data = sel_dat2)

#Data: TOP2A in 9 controls (group asymptomatic) > 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.6898
ci(roc.list$TOP2A)
#95% CI: 0.4228-0.9568 (DeLong)

####################CCNB2###############################
roc.list$CCNB2
#Call:
#  roc.formula(formula = group ~ CCNB2, data = sel_dat2)

#Data: CCNB2 in 9 controls (group asymptomatic) < 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.5417
ci(roc.list$CCNB2)
#95% CI: 0.2708-0.8126 (DeLong)

####################TPX2###############################
roc.list$TPX2
#Call:
#  roc.formula(formula = group ~ TPX2, data = sel_dat2)

#Data: TPX2 in 9 controls (group asymptomatic) < 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.5833
ci(roc.list$TPX2)
#95% CI: 0.3224-0.8442 (DeLong)

####################BUB1###############################
roc.list$BUB1
#Call:
#  roc.formula(formula = group ~ BUB1, data = sel_dat2)

#Data: BUB1 in 9 controls (group asymptomatic) > 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.7269
ci(roc.list$BUB1)
#95% CI: 0.484-0.9697 (DeLong)

####################CDK1###############################
roc.list$CDK1
#Call:
#  roc.formula(formula = group ~ CDK1, data = sel_dat2)

#Data: CDK1 in 9 controls (group asymptomatic) < 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.4815
ci(roc.list$CDK1)
#95% CI: 0.2296-0.7334 (DeLong)

####################KIF2C###############################
roc.list$KIF2C
#Call:
#  roc.formula(formula = group ~ KIF2C, data = sel_dat2)

#Data: KIF2C in 9 controls (group asymptomatic) > 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.6852
ci(roc.list$KIF2C)
#95% CI: 0.4239-0.9465 (DeLong)

####################RRM2###############################
roc.list$RRM2
#Call:
#  roc.formula(formula = group ~ RRM2, data = sel_dat2)

#Data: RRM2 in 9 controls (group asymptomatic) > 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.5556
ci(roc.list$RRM2)
#95% CI: 0.2941-0.817 (DeLong)

####################AURKB###############################
roc.list$AURKB
#Call:
#  roc.formula(formula = group ~ AURKB, data = sel_dat2)

#Data: AURKB in 9 controls (group asymptomatic) < 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.662
ci(roc.list$AURKB)
#95% CI: 0.4132-0.9109 (DeLong)

####################CDCA8###############################
roc.list$CDCA8
#Call:
#  roc.formula(formula = group ~ CDCA8, data = sel_dat2)

#Data: CDCA8 in 9 controls (group asymptomatic) < 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.7685
ci(roc.list$CDCA8)
#95% CI: 0.5552-0.9819 (DeLong)

####################KIF11###############################
roc.list$KIF11
#Call:
#  roc.formula(formula = group ~ KIF11, data = sel_dat2)

#Data: KIF11 in 9 controls (group asymptomatic) > 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.8009
ci(roc.list$KIF11)
#95% CI: 0.5874-1 (DeLong)

####################CCNB1###############################
roc.list$CCNB1
#Call:
#  roc.formula(formula = group ~ CCNB1, data = sel_dat2)

#Data: CCNB1 in 9 controls (group asymptomatic) < 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.5093
ci(roc.list$CCNB1)
#95% CI: 0.2331-0.7855 (DeLong)

####################CDC20###############################
roc.list$CDC20
#Call:
#  roc.formula(formula = group ~ CDC20, data = sel_dat2)

#Data: CDC20 in 9 controls (group asymptomatic) > 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.6343
ci(roc.list$CDC20)
#95% CI: 0.3713-0.8973 (DeLong)

####################UBE2C###############################
roc.list$UBE2C
#Call:
#  roc.formula(formula = group ~ UBE2C, data = sel_dat2)

#Data: UBE2C in 9 controls (group asymptomatic) > 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.5278
ci(roc.list$UBE2C)
#95% CI: 0.2393-0.8163 (DeLong)

####################KIF20A###############################
roc.list$KIF20A
#Call:
#  roc.formula(formula = group ~ KIF20A, data = sel_dat2)

#Data: KIF20A in 9 controls (group asymptomatic) < 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.6111
ci(roc.list$KIF20A)
#95% CI: 0.3028-0.9194 (DeLong)

####################PTTG1###############################
roc.list$PTTG1
#Call:
#  roc.formula(formula = group ~ PTTG1, data = sel_dat2)

#Data: PTTG1 in 9 controls (group asymptomatic) < 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.8796
ci(roc.list$PTTG1)
#95% CI: 0.7326-1 (DeLong)

####################TTK###############################
roc.list$TTK
#Call:
#  roc.formula(formula = group ~ TTK, data = sel_dat2)

#Data: TTK in 9 controls (group asymptomatic) > 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.625
ci(roc.list$TTK)
#95% CI: 0.368-0.882 (DeLong)

####################DLGAP5###############################
roc.list$DLGAP5
#Call:
#  roc.formula(formula = group ~ DLGAP5, data = sel_dat2)

#Data: DLGAP5 in 9 controls (group asymptomatic) > 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.625
ci(roc.list$DLGAP5)
#95% CI: 0.3829-0.8671 (DeLong)

####################BUB1B###############################
roc.list$BUB1B
#Call:
#  roc.formula(formula = group ~ BUB1B, data = sel_dat2)

#Data: BUB1B in 9 controls (group asymptomatic) > 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.6667
ci(roc.list$BUB1B)
#95% CI: 0.3815-0.9519 (DeLong)

####################PBK###############################
roc.list$PBK
#Call:
#  roc.formula(formula = group ~ PBK, data = sel_dat2)

#Data: PBK in 9 controls (group asymptomatic) < 12 cases (group re_detectable_positive_patients).
#Area under the curve: 0.5787
ci(roc.list$PBK)
#95% CI: 0.3615-0.7959 (DeLong)

roc.list <- roc(group ~ CDCA8+KIF11+PTTG1, data = sel_dat2)

g.list <- ggroc(roc.list)
g.list





# with additional aesthetics:
g3 <- ggroc(roc.list, size = 1.2,alpha=.6)
pdf("ROC.pdf", width = 5, height = 4)
g3+ggsci::scale_color_lancet()+theme_bw()
dev.off()

