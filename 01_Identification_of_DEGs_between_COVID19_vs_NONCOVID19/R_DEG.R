setwd("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/01_DEGs")
remove(list=ls())

counts_GSE157103<-read.table("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/DataCollection/COVID19_bulkRNA/GSE157103/GSE157103_genes.ec.tsv",
                      row.names = 1, header = T, sep = "\t")
dim(counts_GSE157103)
#[1] 19472   126

counts_GSE152641<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/DataCollection/COVID19_bulkRNA/GSE152641/GSE152641_Inflammatix_COVID19_counts_entrez.csv",
                      row.names = 1, header = T)
dim(counts_GSE152641)
#[1] 20460    86

library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DESeq2)
counts_GSE152641$geneid <- rownames(counts_GSE152641)
# transform id  
map_dt <- bitr(counts_GSE152641$geneid, fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
dt_merge <- merge(map_dt,counts_GSE152641, by.y = "geneid", by.x = "ENTREZID")
rownames(dt_merge) <- dt_merge$SYMBOL
dt_merge$SYMBOL <- NULL
dt_merge$ENTREZID <- NULL
dim(dt_merge)
#[1] 20437    86
counts_GSE152641<-cbind(dt_merge[,19:80],dt_merge[,1:18],dt_merge[,81:86])
write.csv(counts_GSE152641, file = "/Users/liuzhe/Desktop/cityu/COVID/DataCollection/COVID19_bulkRNA/GSE152641/GSE152641_readcounts.csv", 
          quote = F)
#前62列是COVID19 后24列是healthy controls

counts_GSE152641$geneid<-rownames(counts_GSE152641)
counts_GSE157103$geneid<-rownames(counts_GSE157103)
counts_merged<-merge(counts_GSE152641, counts_GSE157103, by = "geneid")
rownames(counts_merged)<-counts_merged$geneid
counts_merged$geneid<-NULL
write.csv(counts_merged, "discoverydataset_GSE157103_GSE152641.csv", quote = F)
#62列是COVID19 24列是healthy controls 100列是COVID9 26列是healthy controls
dim(counts_merged)
#[1] 19184   212

condition <- c(rep("COVID",62),rep("NONCOVID",24),rep("COVID",100),rep("NONCOVID",26))
batch <- c(rep(1,86),rep(2,126))
sampleTable <- data.frame(condition=condition,batch=batch)
#样本表达数据框列名需要与sampleTable一致。
row.names(sampleTable) <- colnames(counts_merged) 
data = apply(counts_merged, 2, as.integer) ## DESeq2分析需要是整数
row.names(data) <- row.names(counts_merged)

dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleTable, design = ~batch+condition)
dds$condition<- relevel(dds$condition, ref = "NONCOVID") 
dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- vst(dds)
# 这里我还提取了标准化后的表达矩阵，可以用于后续的热图绘制等等
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 
## 4.4定义差异基因
nrDEG <- nrDEG_DESeq2
nrDEG$Group = "notsignificant"
logFC_cutoff <- 1
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange > logFC_cutoff) )] = "upregulated"
nrDEG$Group[which( (nrDEG$padj < 0.05) & (nrDEG$log2FoldChange < -logFC_cutoff) )] = "downregulated"
table(nrDEG$Group)
#downregulated notsignificant    upregulated 
#         125          18123           936 
write.csv(nrDEG,"DEGs.csv",quote = F)


library("ggplot2")
deg<-nrDEG
head(deg)
deg$color<-ifelse(deg$padj<0.05 & abs(deg$log2FoldChange)>logFC_cutoff, ifelse(deg$log2FoldChange< -logFC_cutoff, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
deg$symbol<-row.names(deg)
p <- ggplot(data = deg, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw() +
  xlim(-5,5)
for_label <- deg %>% 
  filter(abs(log2FoldChange) >2& padj < 0.01)
pdf("VP_DEGs_COVID_NONCOVID.pdf")
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )+ theme(legend.position="none")
dev.off()
#riskgenes = c("gene1", "gene2", "gene3", "gene4")
#for_label<-deg[riskgenes,]
normal_gset<-as.data.frame(normal_gset)
deg_down<-subset(deg, Group == "downregulated")
deg_up<-subset(deg, Group == "upregulated")
deg_select<-rbind(deg_down,deg_up)
normalized_exp<-normal_gset[rownames(deg_select),]
exp_deg<-cbind(deg_select,normalized_exp)
write.csv(exp_deg,"deg_exp.csv",quote = F)
exp_deg<-exp_deg[,c(2,5,6,7,10:221)]
dim(exp_deg)
#[1] 1061  216
write.csv(exp_deg,"Table1_COVID_NONCOVID.csv",quote = F)
#exp_deg<-read.csv("Table1_COVID_NONCOVID.csv",header=T, row.names = 1)

rt<-exp_deg[,c(5:66,91:190,67:90,191:216)]
library("pheatmap")

annotation_col = data.frame(Disease = c(rep("COVID19",162),rep("NONCOVID19",50)))
rownames(annotation_col) = colnames(rt)

annotation_row = matrix(exp_deg$Group,nrow=length(exp_deg$Group), ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_deg)
colnames(annotation_row)<-"gene_type"
annotation_row<-as.data.frame(annotation_row)
ann_colors = list(Disease = c(COVID19 = "red", NONCOVID19 = "blue"), 
                  gene_type = c(downregulated = "blue", upregulated = "red"))
#rt1<-t(rt[,2:5])
# only show the top 20 up- and 20 down-regulated genes
heatmap_input<-rbind(rt[1:20,],rt[126:145,])
pdf(file="heatmap_COVID19-NONCOVID19.pdf",width = 6,height = 8)
pheatmap(heatmap_input,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col, 
         annotation_row = annotation_row, gaps_row = 20,gaps_col = 162,angle_col = "45", 
         annotation_colors = ann_colors, show_rownames = T, show_colnames = F, main = "Title")
#fontsize = 1
dev.off()

library("DOSE")
library("org.Hs.eg.db")
library("topGO")
library("clusterProfiler")
library("pathview")
###up DEG functional enrichment analysis
deg_up<-subset(exp_deg, Group == "upregulated")
up_genes <- rownames(deg_up)
test = bitr(up_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
write.csv(summary(go_BP),"GOBP_enrich_up_COVID19-NONCOVID19.csv",row.names =FALSE)
#结果可视化
GOBP_up<-summary(go_BP)
GOBP_up<-GOBP_up[1:30,]
GOBP_up$yvalue<-GOBP_up$Count/850
GOBP_up$xvalue<-GOBP_up$Description
library(ggpubr)
pdf("GOBP_up_COVID19-NONCOVID19.barplot.pdf",height=6,width = 12)
ggbarplot(GOBP_up, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("red","pink"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()
###down DEG functional enrichment analysis
deg_down<-subset(exp_deg, Group == "downregulated")
down_genes <- rownames(deg_down)
test = bitr(down_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  keyType = 'ENTREZID')
write.csv(summary(go_BP),"GOBP_enrich_down_COVID19-NONCOVID19.csv",row.names =FALSE)
#结果可视化
GOBP_down<-summary(go_BP)
GOBP_down<-GOBP_down[1:30,]
GOBP_down$yvalue<-GOBP_down$Count/111
GOBP_down$xvalue<-GOBP_down$Description
pdf("GOBP_down.barplot_COVID19-NONCOVID19.pdf",height=6,width = 12)
ggbarplot(GOBP_down, x = "xvalue", y = "yvalue", orientation = "horiz",
          xlab = "", ylab="Gene Ratio",
          fill = "pvalue",
          color = "white",
          width = 0.6,
          position = position_dodge())+gradient_fill(c("blue","purple"))+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(legend.position = "right")
dev.off()
