rm(list=ls())
## my own data
library(safe)
setwd("F:\\bioinformatics\\module_rank_20191209\\Pathway_selected_fromMsigDB")
data <- as.matrix(read.csv("used_normal_genes.csv",header=T))

used_data <- data[,-c(1,2)]
#genenames <- data[,2]
genenames <- apply(as.matrix(data[,1]),1,function(x)gsub(" ","",x))

##
label <- as.matrix(read.csv("gene_data_lable.csv",header=T))
data_EX <- as.matrix(apply(as.matrix(cbind(used_data[,which(label==1)],used_data[,which(label==0)])),2,function(x)as.numeric(x)))
dim(data_EX)
#rownames(data_EX) <- data[,1]

y <- c(rep(1,length(which(label==1))),rep(0,length(which(label==0))))

#genesets <- readRDS("used_msigdb_pathway_without_Net.RDS")
genesets <- readRDS("small_cluster.RDS")

names(genesets) <- paste("sets",as.character(1:length(genesets)),sep="")

C.mat2 <- getCmatrix(keyword.list = genesets,
                     present.genes = genenames)

results1 <- safe(data_EX, y, C.mat2, print.it = FALSE)
#dd <- safe.toptable(results1, number = 140, description = FALSE)
data <- cbind(as.matrix(names(genesets)),
                    as.matrix(results1@global.pval))
colnames(data) <- c("Setsname","Global")
data <- data[order(as.numeric(data[,2]),decreasing = F),]
#saveRDS(data,"SAFE_Regnet_cluster.RDS")
#write.csv(data,"SAFE_Regnet_cluster.csv",row.names = F,quote = F)

top_10 <- data[1:10,1]
#saveRDS(top_10,"SAFE_top_10_net.RDS")
#top_10 <- readRDS("SAFE_Msig_3.RDS")[1:10,1]
###
nodes <- NULL
for(i in 1:10)
{
  k1 <- as.numeric(gsub("sets","",top_10[i],fixed = T))
  nod1 <- as.matrix(unlist(genesets[[k1]]))
  nodes<- c(nodes,nod1)
}
#352length(unique(nodes))
##986
nodes <- unique(nodes)
library(org.Hs.eg.db)
SY <- mapIds(org.Hs.eg.db,keys=nodes, column="SYMBOL", keytype="ENTREZID", multiVals="first")
head(SY)
length(SY)
write.csv(cbind(nodes,SY),"SAFE_top_10_nodes_RegNet_cluster.csv",
          row.names = F,
          quote = F)
#write.csv(top10,"SAFE_top10.csv",row.names = F,quote = F)
library(infotheo)
used_RNA <- discretize(t(used_data))
MI_cal <- rbind()
## ¼ÆËãMI
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
saveRDS(MI_cal,"SAFE_top_10_MI_Regnet_cluster.RDS")
pdf("SAFE_top_10_MI_boxplot.pdf")
boxplot(MI_cal,col="green")
dev.off()