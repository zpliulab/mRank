rm(list=ls())
setwd("F:\\bioinformatics\\module_rank_20191209\\Pathway_selected_fromMsigDB")
p1 <- "./GSEA/GSEA_check.R"
source(p1)
## GSEA main()
data <- as.matrix(read.csv("used_normal_genes.csv",header=T))

used_data <- data[,-c(1,2)]
genenames <- apply(as.matrix(data[,1]),1,function(x)gsub(" ","",x))
data_used <- as.matrix(apply(as.matrix(data[,-c(1,2)]),2,function(x) as.numeric(x)))
phenotype <- (as.matrix(read.csv("gene_data_lable.csv",header=T)))
rownames(data_used) <- genenames
colnames(data_used) <- colnames(data)[-c(1,2)]
## gene sets
#genesets <- readRDS("used_msigdb_pathway_without_Net.RDS")#readLines("Reg_gene_set_SY_new2.txt")
genesets <- readRDS("small_cluster.RDS")


genesets_sy <- vector("list",length(genesets))
for(i in 1:length(genesets))
{
  nn <- paste("sets",as.character(i),sep="")
  genesets_sy[[i]] <-c(nn,"Non_description",
                       as.matrix(unlist(genesets[[i]])))
}
sep=" "
p2 <- "./GSEA/GSEA_shx1.R"
source(p2)
## 设置初始值
output.directory = ""
doc.string = "GSEA.analysis"
non.interactive.run = F
reshuffling.type = "sample.labels"
nperm = 1000
weighted.score.type = 1
nom.p.val.threshold = -1
fwer.p.val.threshold = -1
fdr.q.val.threshold = 0.25
topgs = 10
adjust.FDR.q.val = F
gs.size.threshold.min = 2
gs.size.threshold.max = 1201
reverse.sign = F
preproc.type = 0
random.seed = 123456
perm.type = 0
fraction = 1.0
replace = F
save.intermediate.results = F
OLD.GSEA = F
use.fast.enrichment.routine = T
kk <- GSEA_main(data_used,genesets_sy,sep)
dim(kk)
saveRDS(kk,"GSEA_Regnet_cluster.RDS")
#write.table(kk,paste(p3,"GSEA_Regnet_cluster.txt",sep=""),sep="\t",row.names = F,quote = F)
## 
#kk <- readRDS("GSEA_path_3.RDS")
## abs(NES) >= 1
## Norm p <= 0.05
## Fdr <= 0.25
as.numeric(as.character(kk[1,7]))
z1 <- intersect(which(as.numeric(kk[,5])>1),which(as.numeric(as.character(kk[,6]))<0.05))
length(z1)
z2 <- intersect(z1,which(as.numeric(as.character(kk[,6]))<0.25))
length(z2)
data1 <-  kk[z1,]
data <- data1[order(abs(as.numeric(data1[,5])),
                    decreasing = T),]
dim(data)
saveRDS(data,"GSEA_Regnet_cluster_top.RDS")
#data <- kk[order(as.numeric(kk[,13])),]
  #k1 <- data[which(data[,13]==0),]
  #gg <- k1[order(k1[,7]),] 
#top_10 <- as.matrix(gg[1:10,1])
top_10 <- as.matrix(data[1:10,1])
#saveRDS(top_10,"GSEA_Regnet_cluster_top.RDS")
#saveRDS(top_10,"GSEA_top_10_new.RDS")
#top_10 <- readRDS("GSEA_top.RDS")
nodes <- NULL
for(i in 1:10)
{
  k1 <- as.numeric(gsub("sets","",top_10[i],fixed = T))
  nod <- as.matrix(unlist(genesets[[k1]]))
  nodes <- c(nodes,nod)
}
nodes <- unique(nodes)
length(nodes)
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

gS <- mapIds(org.Hs.eg.db, keys=nodes, column="SYMBOL", keytype="ENTREZID", multiVals="first")
head(gS)
length(gS)
write.csv(cbind(nodes,gS),"node_in_top_10_GSEA_Regnet_cluster.csv",row.names = F,quote = F)
library(infotheo)
used_RNA <- discretize(t(used_data))#data_used
MI_cal <- rbind()
label<- phenotype
## 计算MI
for(i in 1:length(top_10))
{  
  #i <- 1 
  k1 <- as.numeric(gsub("sets","",top_10[i],fixed = T))
  nodes <- as.matrix(unlist(genesets[[k1]]))
  mmi <- rbind()
  for(j in 1:length(nodes))
  {
    #j <- 1
    k1 <- which(genenames == nodes[j])
    ## Extract split Tumor and Normal samples
    Neg <- c(used_RNA[which(label== 1),k1])
    Neg_l <- c(label[which(label==1)])
    Pos <- c(used_RNA[which(label== 0),k1])
    Pos_l <- c(label[which(label==0)])
    MI <- mutinformation(c(Neg,Pos),c(Neg_l,Pos_l))
    mmi <- rbind(mmi,MI)
  }
  
  MI_cal <- rbind(MI_cal,round(mean(mmi),3))
}
mean(MI_cal)
sd(MI_cal)
saveRDS(MI_cal,"GSEA_top_10_MI_Regnet_cluster.RDS")
pdf("GSEA_top_10_MI_boxplot.pdf")
boxplot(MI_cal,col="green")
dev.off()