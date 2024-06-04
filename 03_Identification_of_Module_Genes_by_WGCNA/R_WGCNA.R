setwd("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/03_WGCNA")
remove(list=ls())

#ref: https://cloud.tencent.com/developer/article/1498123

#clinical data processing
gse157103_cli<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/DataCollection/COVID19_bulkRNA/GSE157103/clinical.csv", header = T)
gse157103_cli_sel<-subset(gse157103_cli, select = c("Sample", "Disease", "Age", "Gender"))
gse157103_cli_sel$Sample_title<-gse157103_cli_sel$Sample
for(i in c(1:126)){
  if(startsWith(gse157103_cli_sel$Sample[i], "C")){
    gse157103_cli_sel$Sample_title[i]<-gse157103_cli_sel$Sample[i]
  }else{
    gse157103_cli_sel$Sample_title[i]<-gsub("N", "NC", gse157103_cli_sel$Sample[i])
  }
}

for(i in c(1:126)){
  if(gse157103_cli_sel$Disease[i] == "COVID"){
    gse157103_cli_sel$Disease[i]<-"COVID19"
  }else{
    gse157103_cli_sel$Disease[i]<-"Healthy control"
  }
}

for(i in c(1:126)){
  if(gse157103_cli_sel$Gender[i] == "male"){
    gse157103_cli_sel$Gender[i]<-"Male"
  }else if(gse157103_cli_sel$Gender[i] == "female"){
    gse157103_cli_sel$Gender[i]<-"Female"
  }else{
    gse157103_cli_sel$Gender[i]<-"NA"
  }
}

gse157103_cli_sel$Age<-gsub("y", "", gse157103_cli_sel$Age)
gse157103_cli_sel$Age<-gsub("\\:", "NA", gse157103_cli_sel$Age)
head(gse157103_cli_sel)
gse157103_cli_sel<-subset(gse157103_cli_sel, select = c("Sample_title", "Disease", "Age", "Gender"))

gse152641_cli<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/DataCollection/COVID19_bulkRNA/GSE152641/clinical.csv", header = T)
gse152641_cli_sel<-subset(gse152641_cli, select = c("Sample_title", "disease", "age", "sex"))
gse152641_cli_sel$disease<-gsub("disease: ", "", gse152641_cli_sel$disease)
gse152641_cli_sel$age<-gsub("age: ", "", gse152641_cli_sel$age)
gse152641_cli_sel$sex<-gsub("Sex: ", "", gse152641_cli_sel$sex)
colnames(gse152641_cli_sel)<-c("Sample_title", "Disease", "Age", "Gender")
head(gse152641_cli_sel)

allsamples_cli<-rbind(gse157103_cli_sel, gse152641_cli_sel)
write.csv(allsamples_cli, "discovery_clinicaldata_GSE157103_GSE152641.csv", quote = F, row.names = F)

cellscore<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/02_ImmuneInfiltraction/Cell_score.csv", 
                    header = T, row.names = 1)
cellscore_t<-t(cellscore)
cellscore_t<-as.data.frame(cellscore_t)
cellscore_t_plasma<-subset(cellscore_t, select ="Plasma cells")
median(cellscore_t_plasma$`Plasma cells`)
#[1] 0.01986445
cellscore_t_plasma$Sample_title<-row.names(cellscore_t_plasma)
merged_data<-merge(cellscore_t_plasma, allsamples_cli, by = "Sample_title")
colnames(merged_data)<-c("Sample_title", "Plasma_cells", "Disease", "Age", "Gender")
cutoff_pla<-median(merged_data$Plasma_cells)
merged_data$Group<-ifelse(merged_data$Plasma_cells>cutoff_pla,"High","Low")
write.csv(merged_data, "discovery_clinicaldata_withGroup_GSE157103_GSE152641.csv", quote = F, row.names = F)


#expression profile data processing
FPKM<-read.csv("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/02_ImmuneInfiltraction/FPKM_discoverydataset_GSE157103_GSE152641.csv", 
               header = T, row.names = 1)
exp<-FPKM
exp_t<-t(exp)
exp_t<-as.data.frame(exp_t)
exp_t$Sample_title<-row.names(exp_t)

#expression profile with group label
#################################step1_构建WGCNA所需表达矩阵#####################################################################
FPKM_group<-merge(exp_t, merged_data, by = "Sample_title")
row.names(FPKM_group)<-FPKM_group$Sample_title
FPKM_group$Sample_title<-NULL

group_list<-factor(FPKM_group$Group)
table(group_list)
#group_list
#High  Low 
#106  106 
write.csv(FPKM_group, "discovery_FPKM_clinical_withGroup_GSE157103_GSE152641.csv", quote = F)

disease_list<-factor(FPKM_group$Disease)
table(disease_list)
disease_list
#COVID19 Healthy control 
#   162              50

#https://www.jianshu.com/p/25905a905086

datTraits <- data.frame(row.names=rownames(FPKM_group),
                        plasma_cells=group_list,
                        covid_19=disease_list)
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
#enableWGCNAThreads()
disableWGCNAThreads()
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
datExpr <- FPKM_group[1:212, 1:19060]
nSamples = nrow(datExpr)
nSamples
#[1] 212
datExpr_t<-t(datExpr)
dim(datExpr_t)
#[1] 19060   212

###########基因水平的过滤
## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(datExpr_t,1,mad)
dataExprVar <- datExpr_t[which(m.mad > 
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))
dim(dataExpr)
#[1]   212 14295
dim(datExpr_t)
#[1] 19060   212
## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
#[1]   212 14295
head(dataExpr)[,1:8]
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
par(mfrow = c(1,1))
pdf("SampleTree.pdf",width = 20,height = 5)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
abline(h=150000,col="red")
clust = cutreeStatic(sampleTree, cutHeight = 150000, minSize = 10)
table(clust)
#clust
#0   1 
#3 209
keepSamples = (clust==1)
datExpr = dataExpr[keepSamples,]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dataExpr<-datExpr

#################################step2_找到合适的beta值#####################################################################
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
RpowerTable=pickSoftThreshold(datExpr, powerVector=powers)[[2]]
cex1=0.7
pdf(file="softThresholding.pdf")
par(mfrow=c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
text(RpowerTable[,1], RpowerTable[,5], labels=powers, cex=cex1,col="red")
dev.off()
power = sft$powerEstimate
power
#[1] 1
power<-6



# 检验选定的β值下记忆网络是否逼近 scale free

# 获得临近矩阵：
softPower <- 6
adjacency = adjacency(dataExpr, power = softPower);
# 将临近矩阵转为 Tom 矩阵
TOM = TOMsimilarity(adjacency);
# 计算基因之间的相异度
dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average");

# 基因多的时候使用下面的代码：
k <- softConnectivity(dataExpr,power=softPower) 
pdf("scalefree_softpower.pdf")
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()
#可以看出k与p(k)成负相关(相关性系数0.85),说明选择的β值能够建立基因无尺度网络。


##################################step3_对基因聚类成模块并且可视化#############################################################
##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0("COVID19", ".tom"),
                       verbose = 3)
# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。 
table(net$colors)
#  0    1    2    3    4    5    6    7    8    9 
#657 8928 1583 1173  794  472  214  197  146  131

## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
table(moduleColors)
#moduleColors
#black      blue     brown     green      grey   magenta      pink       red turquoise    yellow 
#  197      1583      1173       472       657       131       146       214      8928       794
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
pdf("cluster_dendrograms.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


#######################step4_将样本信息添加进去,分析样本性状与模块的关系########################################
design1=model.matrix(~0+ datTraits$plasma_cells)
design1 = as.data.frame(design1)
colnames(design1)=levels(datTraits$plasma_cells)

design2=model.matrix(~0+ datTraits$covid_19)
design2 = as.data.frame(design2)
colnames(design2)=levels(datTraits$covid_19)

design<-cbind(design1, design2)
rownames(design)<-rownames(datTraits)

MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
pdf("moduleTraitCor.pdf")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


############################step5_感兴趣性状的模块的具体基因分析###############################################
## 首先计算模块与基因的相关性矩阵
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="")

## 再计算性状与基因的相关性矩阵
geneTraitSignificance = as.data.frame(cor(dataExpr, design, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(design), sep="")
names(GSPvalue) = paste("p.GS.", names(design), sep="")


#只保留与Plasma cells相关性大于0.2 并且与COVID19相关性也大于0.2的模块
## 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = c("black")
module
#[1] "black"

column = match(module, modNames);
moduleGenes = moduleColors==module;


##############################step6_网络的可视化#####################################################
# -----------------------------------------------------------------------------------
## 首先针对所有基因画热图
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
geneTree = net$dendrograms[[1]]; 
dissTOM = 1-TOMsimilarityFromExpr(dataExpr, power = 6); 
plotTOM = dissTOM^7; 
diag(plotTOM) = NA; 
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

## 最后画模块和性状的关系
# Recalculate module eigengenes
MEs = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MET = orderMEs(cbind(MEs, design))
# Plot the relationships among the eigengenes and the trait
pdf("Relationship_eigengenesAndtraits.pdf")
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
# Plot the dendrogram
dev.off()
pdf("Relationship_modelAndtraits.pdf")
## 模块与性状的聚类图
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)

pdf("heatmap_modelAndtraits.pdf")
## 性状与模块热图
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

################################step7_提取指定模块的基因名################################################
# --------------------------------------------------------------------------------------
# Select module
#module = c("black", "grey")
#black 197
#grey 657
#197+657=854

# Select module probes
probes = colnames(dataExpr) ## 我们例子里面的probe就是基因名
#inModule = (moduleColors %in% module);
inModule <- (moduleColors %in% c("black"))
modProbes = probes[inModule]; 
head(modProbes)
#[1] "ABCB9"   "ALG14"   "ALG5"    "APOO"    "AURKB"   "B4GALT2"
length(modProbes)
#[1] 197

#######################Step9: 模块的导出##########################################
# --------------------------------------------------------------------------------------
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(dataExpr, power = 6); 
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
)
save(cyt, file = "Cyt.RData")
class(cyt)
#[1] "list"
dim(cyt$edgeData)
#[1] 14369     6
dim(cyt$nodeData)
#[1] 197   3

edg_data<-cyt$edgeData
edg_data<-subset(edg_data, select = c("fromNode", "toNode", "weight"))
head(edg_data)
write.csv(edg_data, "edge_table.csv", quote = F, row.names = F)

nod_data<-cyt$nodeData
nod_data<-subset(nod_data, select = c("nodeName", "nodeAttr[nodesPresent, ]"))
head(nod_data)
colnames(nod_data)<-c("nodeName", "module")
write.csv(nod_data, "node_table.csv", quote = F, row.names = F)


##########Correlation between module membership and gene significance of gene in selected module
#Scatter plots of modules eigengenes in the black module.
module = c("black")
column = match(module, modNames)
modulesGene = moduleColors==module
pdf("Model_gene_corr_p.pdf")
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlabs = paste("Module Membership in", module, "module"),
                   ylabs = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()




