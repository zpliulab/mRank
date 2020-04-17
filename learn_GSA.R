###############
##pathway from KEGG, Biocarta, Reactome,GSA my own data try
###############
rm(list=ls())
library(GSA)
setwd("F:\\bioinformatics\\module_rank_20191209\\Pathway_selected_fromMsigDB")
data <- as.matrix(read.csv("used_normal_genes.csv",header=T))

used_data <- data[,-c(1,2)]

genenames <- apply(as.matrix(data[,1]),1,function(x)gsub(" ","",x))

##
label <- as.matrix(read.csv("gene_data_lable.csv",header=T))
data_EX <- as.matrix(apply(as.matrix(cbind(used_data[,which(label==1)],used_data[,which(label==0)])),2,function(x)as.numeric(x)))
dim(data_EX)

y <- c(rep(1,length(which(label==1))),rep(2,length(which(label==0))))


#genesets <- readRDS("used_msigdb_pathway_without_Net.RDS")
genesets <- readRDS("small_cluster.RDS")
names(genesets) <- paste("sets",as.character(1:length(genesets)),sep="")

## 需要设置最大和最小的geneset 的size
## 5 1201
GSA.obj <- GSA(data_EX,y, genenames=genenames, genesets=genesets,minsize=5,maxsize=1201, resp.type="Two class unpaired", nperms=100)
#mode(GSA.obj)
#setwd("./GSA/learn")
#getwd()
GSA.listsets(GSA.obj, geneset.names=names(genesets),FDRcut=.5)
pdf("GSA_path.pdf")
GSA.plot(GSA.obj)
dev.off()
data <- cbind(names(genesets),abs(as.numeric(GSA.obj[["GSA.scores"]])))
dim(data)
dd <- data[order(as.numeric(data[,2]),decreasing=T),]
write.csv(dd,"GSA_Regnet_cluster.csv",row.names = F,quote = F)
#dd <- as.matrix(read.csv("GSA_path_3.csv",header=T))
top_10 <- dd[1:10,1]
#saveRDS(top_10,"GSA_top_10_Regnet.RDS")
nodes <- NULL
for(i in 1:10)
{
  k1 <- as.numeric(gsub("sets","",top_10[i],fixed = T))
  nod1 <- as.matrix(unlist(genesets[[k1]]))
  nodes <- c(nodes,nod1)
}
nodes <- unique(nodes)
length(nodes)
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

S <- mapIds(org.Hs.eg.db, keys=nodes, column="SYMBOL", keytype="ENTREZID", multiVals="first")
head(S)
length(S)
write.csv(cbind(nodes,S),"nodes_in_top_10_GSA_Regnet_cluster.csv",row.names = F,quote = F)
library(infotheo)
used_RNA <- discretize(t(used_data))
dim(used_RNA)
MI_cal <- rbind()
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
saveRDS(MI_cal,"GSA_top_10_MI_Regnet_cluster.RDS")
pdf("GSA_top_10_MI_boxplot.pdf")
boxplot(MI_cal,col="red")
dev.off()
## 计算AUC


#da <- as.matrix(GSA.obj[["GSA.scores"]])
#length(da)
#capture.output(da, file = "GSA_store_results.txt")
#oo <- GSA.listsets(GSA.obj,geneset.names=names1,FDRcut = 1)

#zz <- cbind(names1,as.numeric(da))
#dim(zz)
#colnames(zz) <- c("gene_sets_name","Score")
#write.csv(cbind(names1,da), file = "GSA_store_results.csv",quote=F,row.names=F)
#k1 <- oo[["negative"]]
#k2 <- oo[["positive"]]
#write.csv(k1, file = "GSA_store_results_ne.csv",quote=F,row.names=F)
#write.csv(k2, file = "GSA_store_results_po.csv",quote=F,row.names=F)

#capture.output(oo[["negative"]], file = "GSA_store_results_ne.txt",quote=F,row.names=F)
#capture.output(oo[["positive"]], file = "GSA_store_results_po.txt",quote=F,row.names=F)
#capture.output(rbind(oo[["negative"]],oo[["positive"]]), file = "GSA_store_results_Final.txt",quote=F,row.names=F)
