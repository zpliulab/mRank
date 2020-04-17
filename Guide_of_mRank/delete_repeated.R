rm(list=ls())

setwd("F:\\bioinformatics\\module_rank_20191209\\20200319\\Guide_of_mRank")
p1 <- "./Module_store/"
file  <- list.files(p1)
length(file)
p2 <- "F:\\bioinformatics\\module_rank_20191209\\20200319\\Guide_of_mRank\\"
DEGs <- as.numeric(as.matrix(read.table(paste(p2,"DE_genes/DE_gene_with_ID.txt",sep=""),header=F,sep="\t"))[,3])

DE_AUC <- rbind()

for(i in 1:length(file))
{
     #i <- 300
    kk <- file[i]
    k1 <- as.matrix(unlist(strsplit(kk,".",fixed = T)))
   k11 <- k1[2]
  name <- as.matrix(unlist(strsplit(k1[1],"_",fixed = T)))[1]
    if(k11 == "RDS")
    {
        data <- readRDS(paste(p1,kk,sep=""))
        all_gene <- as.numeric(data)
        kk1 <- intersect(data,as.character(DEGs))
        #length(data)
      DE_AUC <- rbind(DE_AUC,c(name,length(kk1),(length(data)-1),data[length(data)-1],as.numeric(data[length(data)])))
    }
}
#saveRDS(DE_AUC,"2608_info.RDS")
### delete repeated modules
### if one module have the same DEG gene sets with another module 
### we choose the module with larger AUC值
Lable <- NULL
AUC <- NULL
count <- 1
for(i in 1:(length(file)-1))
{
  #i  <- 1
   kk <- file[i]
   k1 <- as.matrix(unlist(strsplit(kk,".",fixed = T)))
  k11 <- k1[2]
 name <- as.matrix(unlist(strsplit(k1[1],"_",fixed = T)))[1]
  if(k11 == "RDS")
  {
    data <- readRDS(paste(p1,kk,sep=""))
    gene <- data[-(length(data))]
      zz <- cbind()
     auc <- cbind()
     for(jj in (i+1):length(file))
     {
       kk2 <- file[jj]
        k2 <- as.matrix(unlist(strsplit(kk2,".",fixed=T)))
       k22 <- k2[2]
       name1 <- as.matrix(unlist(strsplit(k2[1],"_",fixed=T)))[1]
       if(k22 == "RDS")
       {
         data2 <- readRDS(paste(p1,kk2,sep=""))
         gene2 <- data2[-length(data2)]
         inn <- intersect(gene,gene2)
         if(length(inn)==length(gene))
         {
           zz <- cbind(zz,name1)
          auc <- cbind(auc,data2[length(data2)] )
         }
       }
     }
      if(length(zz)>0)
        {
        Lable[[count]] <- c(name,zz)
          AUC[[count]] <- c(data[length(data)],auc)
                 count <- count+1
       }
  }
}

#saveRDS(Lable,"repeat_DE.RDS")
#saveRDS(AUC,"repeat_DE_AUC.RDS")

##筛选重复里AUC最大的
DE_in <- NULL
DE_final_chosen <- rbind()
for(i in 1:length(AUC))
{
  #i <- 1
  k1 <- as.matrix(unlist(AUC[[i]]))
  k2 <- as.matrix(unlist(Lable[[i]]))
  used_k1 <- which(as.numeric(k1)==max(as.numeric(k1)))[1]
  DE_in <- c(DE_in,k2)
  DE_final_chosen <- rbind(DE_final_chosen,k2[used_k1])
}
#saveRDS(DE_in,"repeat_DE_Label.RDS")
#saveRDS(DE_final_chosen,"repeat_final_chose.RDS")
##先取差集再取并集
all_used_DE <-rbind()
for(j in 1:length(file))
{
  #j  <- 1
  kk <- file[j]
  k1 <- as.matrix(unlist(strsplit(kk,".",fixed = T)))
  #k11 <- k1[2]
  name <- as.matrix(unlist(strsplit(k1[1],"_",fixed = T)))[1]
  all_used_DE <-rbind(all_used_DE,name)
}

##
DE_list <- union(setdiff(all_used_DE,DE_in),DE_final_chosen)
length(DE_list)
##2561
saveRDS(DE_list,"DE_list_final.RDS")
### 
library(pracma)

Final_module <- matrix(0,nrow=2561,ncol=20)
count <- 1
for(i in 1:length(file))
{
  kk <- file[i]
  k1 <- as.matrix(unlist(strsplit(kk,".",fixed = T)))
  k11 <- k1[2]
  name <- as.matrix(unlist(strsplit(k1[1],"_",fixed = T)))[1]
  z1 <- which(name %in% DE_list)
  if(k11 == "RDS")
  {
    if(isempty(z1)==F)
    {
     data <- readRDS(paste(p1,kk,sep=""))
     all_gene <- as.matrix(data)
     kk1 <- intersect(data,as.character(DEGs))
     #length(data)
     Final_module[count,1:length(all_gene)] <- all_gene
     count <- count+1
    }
  }
}
#saveRDS(Final_module,"module_final.RDS")
write.csv(Final_module,"module_final.csv",quote = F,row.names = F)
