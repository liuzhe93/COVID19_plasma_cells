setwd("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/06_PPI")
remove(list=ls())

#STRING: 11.5 online version
net_data<-read.table("string_interactions_short.tsv", header = T)
dim(net_data)
#[1] 2443   13
net_data_fil<-net_data[which(net_data$combined_score>0.8),]
dim(net_data_fil)
#[1] 1283  13
length(unique(c(net_data_fil$node1, net_data_fil$node2)))
#[1] 89
#the number of edges: 1283
#the number of nodes: 89
write.csv(net_data_fil, row.names = F, "network.csv", quote = F)



#cytoscape: 3.9.0
top20<-c("CCNA2", "TOP2A", "CCNB2", "TPX2", "BUB1", "CDK1", "KIF2C", "RRM2", "AURKB", "CDCA8", 
         "KIF11", "CCNB1", "CDC20", "UBE2C", "KIF20A", "PTTG1", "TTK", "DLGAP5", "BUB1B", "PBK")
degs<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/01_DEGs/DEGs.csv", header = T, row.names = 1)
degs_sel<-degs[top20,]
degs_sel
degs_exp<-read.csv("c:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/01_DEGs/deg_exp.csv", header = T, row.names = 1)
exp_sel<-degs_exp[top20,]
exp_sel
row_mean<-apply(exp_sel[,10:221], 1, mean)
row_mean<-as.data.frame(row_mean)
exp_sel$AverageExp<-mean(exp_sel[,10:221])
data_select<-subset(exp_sel, select = c("log2FoldChange", "pvalue", "padj", "Group"))
mydata<-cbind(data_select, row_mean)
write.csv(mydata, "top20_detail.csv", quote = F)

