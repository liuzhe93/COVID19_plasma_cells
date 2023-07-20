setwd("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/05_SharedGenes")
remove(list=ls())

deg<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/01_DEGs/DEGs.csv")
wgcna<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/03_WGCNA/node_table.csv")
scRNA<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/04_scRNA/covid19_subset_plasma_cells.markers.csv")


deg_sel<-deg[which(deg$Group!="notsignificant"),]
genes_deg<-deg_sel$X
length(genes_deg)
#[1] 1061
genes_wgcna<-wgcna$nodeName
length(genes_wgcna)
#[1] 197
genes_scRNA<-rownames(scRNA)
length(genes_scRNA)
#[1] 296

overlap_deg_wgcna<-intersect(genes_deg, genes_wgcna)
overlap_deg_scRNA<-intersect(genes_deg, genes_scRNA)
overlap_wgcna_scRNA<-intersect(genes_wgcna, genes_scRNA)
overlap_deg_wgcna_scRNA<-intersect(overlap_deg_wgcna, genes_scRNA)

length(overlap_deg_wgcna)
#[1] 108
length(overlap_deg_scRNA)
#[1] 8
length(overlap_wgcna_scRNA)
#[1] 16
length(overlap_deg_wgcna_scRNA)
#[1] 5
dat<-c("DEGs" = 950, 
       "Module_genes" = 78, 
       "Marker_genes" = 277,
       "DEGs&Module_genes" = 103, 
       "DEGs&Marker_genes" = 3, 
       "Module_genes&Marker_genes" = 11,
       "DEGs&Module_genes&Marker_genes" = 5)

library("eulerr")
pdf("shared_genes.pdf")
plot(euler(dat),
     fills = list(fill=c("red","blue",
                         "green","darkgreen",
                         "orange","black",
                         "purple"),
                  alpha=0.5),
     quantities = list(c(950,78,277,
                    103,3,11,5),
                    col="black",
                    cex=4),
     labels = list(col="white",font=3,cex=2),
     edges = list(col="darkgreen",lwd=5,
                  lty=1:3), main = list(label=c("Shard_genes"),cex=5),
     legend = list(labels=c("DEGs","Module_genes","Marker_genes"),
                   cex=1))
dev.off()

library("UpSetR")

upset_list<-list(genes_deg, genes_wgcna, genes_scRNA)
names(upset_list)<-c("DEGs", "Module_genes", "Marker_genes")
pdf("Upset_shared.pdf", width = 9, height = 6)
upset(fromList(upset_list),
      nsets = 100,
      nintersects = 40,
      order.by = "freq",
      keep.order = F,
      mb.ratio = c(0.6, 0.4),
      text.scale = 2)
dev.off()

twice_genes<-c(genes_deg, genes_wgcna, genes_scRNA)
genes_counts<-table(twice_genes)
genes_data<-as.data.frame(genes_counts)
genes_data_fil<-genes_data[which(genes_data$Freq>=2),]
shared_genes<-genes_data_fil$twice_genes
write.csv(shared_genes, "shared_genes.csv", quote = F, row.names = F)


#reference
#https://blog.csdn.net/weixin_42655515/article/details/113830462
#https://www.cnblogs.com/yanjiamin/p/12122215.html


library("ggpubr")
library("DOSE")
library("org.Hs.eg.db")
library("topGO")
library("clusterProfiler")
library("pathview")
test = bitr(shared_genes, #dataset
            fromType="SYMBOL",
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db") 
head(test,2)
GO<-enrichGO(test$ENTREZID,
             keyType = 'ENTREZID',
             OrgDb = org.Hs.eg.db,
             ont="ALL",
             pAdjustMethod = 'BH',
             pvalueCutoff = 1, 
             qvalueCutoff = 1,
             readable = T)
go<-as.data.frame(GO)
View(go)
table(go[,1]) #查看BP,CC,MF的统计数目
#  BP   CC   MF 
#1841  157  263
write.csv(go, "GOBP_enrich_readable.csv")

go_MF<-go[go$ONTOLOGY=="MF",][1:10,]
go_CC<-go[go$ONTOLOGY=="CC",][1:10,]
go_BP<-go[go$ONTOLOGY=="BP",][1:10,]
go_enrich_df<-data.frame(ID=c(go_BP$ID, go_CC$ID, go_MF$ID),
                         Description=c(go_BP$Description, go_CC$Description, go_MF$Description),
                         GeneNumber=c(go_BP$Count, go_CC$Count, go_MF$Count),
                         type=factor(c(rep("biological process", 10), rep("cellular component", 10),rep("molecular function",10))),
                         levels=c("molecular function", "cellular component", "biological process"))

library("ggrepel") # 标签相关
library("stringr") # 标签换行

go_enrich_df$Description <- factor(go_enrich_df$Description,levels = rev(go_enrich_df$Description))

go_bar <- ggplot(data = go_enrich_df, # 绘图使用的数据
                 aes(x = Description, y = GeneNumber,fill = type))+ # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
ggsave(go_bar,filename = "GO_Barplot.pdf",width = 9,height = 7)

kegg_data<-read.table("KEGG_fromDAVID/david.ncifcrf.gov_data_download_chart_C7F4CE1FE3D71689320839779.txt", header = T, sep = "\t")
kegg_data_fil<-kegg_data[which(kegg_data$PValue<0.05),]
kegg_data_srt<-kegg_data_fil[order(kegg_data_fil$PValue),]

kegg_data_srt$Term <- factor(kegg_data_srt$Term,levels = rev(kegg_data_srt$Term))



kegg_bar <- ggplot(data = kegg_data_srt, # 绘图使用的数据
                 aes(x = Term, y = Count, fill = PValue))+ # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "KEGG terms",y = "GeneNumber",title = "Barplot of Enriched KEGG Terms")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距

ggsave(kegg_bar,filename = "KEGG_Barplot.pdf",width = 9,height = 4)

