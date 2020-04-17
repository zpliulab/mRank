###
rm(list=ls())
setwd("F:/bioinformatics/module_rank_20191209/Crank_preapre")
#Net input
Net <- as.matrix(read.table("renewed_RegNetwork.txt",sep="\t",header = T))
#dim(Net)
## Gene EX 
data <- as.matrix(read.csv("used_normal_genes.csv",header=T))

used_data <- data[,-c(1,2)]
#genenames <- data[,2]
genenames <- apply(as.matrix(data[,1]),1,function(x)gsub(" ","",x))
label <- as.matrix(read.csv("gene_data_lable.csv",header=T))
data_EX <- as.matrix(apply(as.matrix(cbind(used_data[,which(label==1)],used_data[,which(label==0)])),2,function(x)as.numeric(x)))
dim(data_EX)

k1 <- which(Net[,1] %in% as.numeric(genenames))
k2 <- which(Net[,2] %in% as.numeric(genenames))
kk <- intersect(k1,k2)
length(kk)#120647
#saveRDS(Net[kk,],"used_RegNet.RDS")
rm(list=ls())
setwd("F:/bioinformatics/module_rank_20191209/Crank_preapre")

Net <- readRDS("used_RegNet.RDS")
library(igraph)

PP <- graph.data.frame(Net,directed = F)
sim_PP <- igraph::simplify(PP,remove.loops = T,remove.multiple = T)
fg <- cluster_fast_greedy(sim_PP)
sizes(fg)
Com <- function(x)
{
  library(igraph)
  k11 <- which(Net[,1] %in% as.numeric(x))
  k21 <- which(Net[,2] %in% as.numeric(x))
  kk1 <- intersect(k11,k21) 
  PP1 <- graph.data.frame(Net[kk1,],directed = F)
  sim_PP1 <- igraph::simplify(PP1,remove.loops = T,remove.multiple = T)
  fg1 <- cluster_fast_greedy(sim_PP1)
  #kk1 <- membership(fg1)
  #sizes(fg1)
  return(fg1)
}
Com_sp <- function(fg)
{
     Coom <- NULL
       c1 <- 1
       ST <- NULL
       c2 <- 1
    for(i in 1:length(fg))
     {
       if(length(fg[[i]])>20)
        {
           k1 <- Com(fg[[i]])
           for(j in 1:length(k1))
           {
             if(length(k1[[j]])>20)
              {
               Coom[[c1]] <- as.matrix(unlist(k1[[j]]))
               c1 <- 1+c1
             }
               else{
                 ST[[c2]] <- as.matrix(unlist(k1[[j]]))
                 c2 <- 1+c2
               }
           }
       }else{
         ST[[c2]] <- as.matrix(unlist(fg[[i]]))
               c2 <- 1+c2
       }
    }

 return(Coom)
}
Com_sp1 <- function(fg)
{
  Coom <- NULL
  c1 <- 1
  ST <- NULL
  c2 <- 1
  for(i in 1:length(fg))
  {
    if(length(fg[[i]])>20)
    {
      k1 <- Com(fg[[i]])
      for(j in 1:length(k1))
      {
        if(length(k1[[j]])>20)
        {
          Coom[[c1]] <- as.matrix(unlist(k1[[j]]))
          c1 <- 1+c1
        }
        else{
          ST[[c2]] <- as.matrix(unlist(k1[[j]]))
          c2 <- 1+c2
        }
      }
    }else{
      ST[[c2]] <- as.matrix(unlist(fg[[i]]))
      c2 <- 1+c2
    }
  }
  
  return(ST)
}
##
YY <- Com_sp(fg)#第一次之后需要
ST <- Com_sp1(fg)#第一次之后留下的
YY1 <- Com_sp(YY)#第二次之后需要
ST1 <- Com_sp1(YY)#第二次之后留下的
YY2 <- Com_sp(YY1)#第三次之后需要的
ST2 <- Com_sp1(YY1)#第三次之后留下的
all_cluster <- c(ST,ST1,ST2,YY2)

All_cluster1 <-NULL
cc <-1
#dir.create("./small_cluster")
for(i in 1:length(all_cluster))
{
  if(length(all_cluster[[i]])>4)
  {
    #All_cluster1[[cc]] <- all_cluster[[i]]
                 
    d1 <- c(paste("sets",as.character(cc),sep=""),
            as.matrix(unlist(all_cluster[[i]])))
    write.table(t(d1),paste(paste("./small_cluster/",as.character(cc),sep=""),
                        ".txt",sep=""),
                 row.names=F,quote=F)
    cc <- 1+cc
  }
}
#saveRDS(All_cluster1,"small_cluster.RDS")
###最后整合一下
##Coom3
## ST1 ST2 ST3
##筛选node >4
#kk <- membership(fg)
### 下面是聚类，
###之后是把大的分成小的
###把网络结构和module 一起输入到CRank中，得到多个方法的
### top 10 的排序
###Cluster by edgebetwenness
#eb <- cluster_edge_betweenness(sim_PP)



if(F){
length(g1)
CM <- vector("list",2500)
for(j in 1:length(g1))
{
  if(length(g1[[j]])<20)
  {
    if(length(g1[[j]])>4)
    {
      
    }
  }
}
z1 <- sizes(g1)

  k11 <- which(Net[,1] %in% as.numeric(z1))
k21 <- which(Net[,2] %in% as.numeric(z1))
kk1 <- intersect(k11,k21) 
PP1 <- graph.data.frame(Net[kk1,],directed = T)
#PP1 <- igraph::simplify(PP1,remove.loops = T,remove.multiple = T)
}
if(F){
ifo <- cluster_infomap(sim_PP)
membership(ifo)
sizes(ifo)
wt <- cluster_walktrap(sim_PP,steps = 4)
membership(wt)
length(unique(sizes(wt)))
lg <- cluster_leading_eigen(sim_PP)
membership(lg)
sizes(lg)
length(unique(sizes(lg)))
cmo <- cluster_optimal(sim_PP)
membership(cmo)
length(unique(sizes(cmo)))
}