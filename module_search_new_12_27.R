

rm(list=ls())
setwd("D:\\shx_bioinformatics\\Module_rank\\module_rank_regnetwork\\19_12_13_new_procedure")

#######################################################
## 实现无放回抽样7次
## store
resample1 <-rbind()
all_loc <- c(1:371)
chosen_already <- NULL

for(bb in 1:7)
  {
   set.seed(1234)
   res<- sample(1:(371-(bb-1)*50),50) 
   all_loc <- setdiff(all_loc,chosen_already)
   chosen <- as.matrix(all_loc[res])
   resample1 <- rbind(resample1,t(as.matrix(chosen)))
   chosen_already <- c(chosen_already,as.numeric(chosen))
  }
#####
## SVM
##### 

Auc_all <- function(Nodes){

  SVM_leave_one_out <-function(mdata,label)
 {
  library("e1071")
  library(pROC)
  svm <- NULL
  for(i in 1:nrow(label))
  {
   mtrain <- label[-i,]#训练数据标签
traindata <- mdata[-i,]#训练数据
  predata <- mdata[i,]#测试数据
    model <- svm(traindata,mtrain,type="C-classification",kernel="radial")#留一法训练SVM 分类
     pred <- predict(model,t(predata),decision.values=TRUE)#预测值
        f <- attr(pred,"decision.values")#提取决策值
      prb <- 1/(1+exp(-f))#把决策值使用sigmoid 函数映射到0~1之间
      svm <- rbind(svm,prb)
  }
  return (svm)
  }

    nodes <- as.matrix(unlist(Nodes))
      loc <- which (G_ID %in% nodes) 
      auc <- cbind()
   for(kk in 1:7)
   {
     De_data_T <- RT_data[resample1[kk,],loc]
     De_data_N <- RN_data[,loc]
    Data_label <- as.matrix(c(rep(1,50),rep(0,50)))
      mdata_DE <- as.matrix(rbind(as.matrix(De_data_T),as.matrix(De_data_N)))
    dim(mdata_DE)
    result1 <- SVM_leave_one_out(mdata_DE,Data_label)
      show1 <- roc(as.factor(Data_label),as.vector(result1),plot=F,print.auc=T,print.thres=T,col="6",lty=1)
         mn <- as.numeric(show1[["auc"]]) 
        auc <- cbind(auc,round(mn,4))
   }
   store <- list(c(nodes,round(sum(auc)/7,4)))
   return(store)
}
###################################################################




##################################################################
#### 网络数据及RNA-seq 数据及其label以及DEGs 
p1 <- "D:/shx_bioinformatics/Module_rank/module_rank_regnetwork/"

GRN <- as.matrix(read.table(paste(p1,"Regnetwork/renewed_RegNetwork_DMI_new.txt",sep=""),header=T,sep=" "))


normal_genes <- as.matrix(read.csv(paste(p1,"RNA_seq/used_normal_genes.csv",sep=""),header=T))
G_ID <- as.numeric(normal_genes[,1])
R_data <- normal_genes[,3:dim(normal_genes)[2]]
R_data <- apply(as.matrix(R_data),1,function(x) as.numeric(x))
Lable <-as.matrix(read.csv(paste(p1,"RNA_seq/gene_data_lable.csv",sep=""),header=T))

TT <- which(Lable==1)#dim(as.matrix(TT))
NN <- which(Lable==0)
RT_data <- R_data[TT,]
RN_data <- R_data[NN,]

DEGs <- as.numeric(as.matrix(read.table(paste(p1,"DE_genes/DE_gene_with_ID.txt",sep=""),header=F,sep="\t"))[,3])
length(DEGs)#2850
#######################
##程序部分（1）
##创建graph图对象
#######################
library(igraph)
REG <- graph.data.frame(GRN[,1:2],directed=F)
Nodes <- as.matrix(get.vertex.attribute(REG)[[1]])
Module_search <- NULL
DE_used <- rbind()
###选取当前 DE 放进去
FFDE <- setdiff(DEGs,DE_used)
AUCCC <- rbind() 
for(mmmm in 1:length(FFDE))
{
DE_cur <- as.character(FFDE[mmmm])

## 查看DEGs中是否有与当前DE的最短路径的DEG

other_DE_Nodes <- setdiff(as.character(DEGs),DE_cur) #length(other_DE_Nodes)
DE_used <- NULL
Path_short <- NULL 
Size_path <- rbind()# path 的length
library(pracma)
count <- 1
other_DEs_with_paths <- rbind()
for(j in 1:length(other_DE_Nodes))
{
   #j <- 1
   if(isempty(which(DE_cur == Nodes))==F)
    {  
            if(isempty(which(other_DE_Nodes[j]== Nodes))==F)
             {
           Path_cur <- shortest_paths(REG, DE_cur, to = other_DE_Nodes[j])
                        if(length(as.matrix(unlist(Path_cur$vpath)))>0)
                            {
                              Path_short[[count]] <- list(Nodes[as.matrix(unlist(Path_cur$vpath))])
                                        Size_path <- rbind(Size_path,length(as.matrix(unlist(Path_cur$vpath))))
                                            count <- count +1
                             other_DEs_with_paths <- rbind(other_DEs_with_paths,other_DE_Nodes[j])

                            }
            }
     }
}
#length(Size_path)
#length(Path_short)
#length(other_DEs_with_paths)
#######
##寻找最短的short path
#######
 Min_Size_path <- min(Size_path[Size_path>0])
 Loc <- which(Size_path == Min_Size_path)#

####### 最短short path的DEGs
other_DE_exist <-  other_DEs_with_paths[Loc]



###############
## part I 
###############
Base_gene <- unique(c(DE_cur,other_DE_exist))
#DE_used <- rbind(DE_used,Base_gene)
### DEG_cur及其他DE的一阶邻居
First_NE <- rbind()

for(kk in 1:length(Base_gene))
{
  First_NE1 <- setdiff(Nodes[as.matrix(unlist(ego(REG, order = 1, Base_gene[kk])))],Base_gene[kk])
  First_NE <- unique(union(First_NE,unique(First_NE1)))
}
#count_First_NE <- length(unique(First_NE))
### SVM RFE 筛选top 50
### 写成子函数

library(caret)
First_NE <- setdiff(First_NE,Base_gene)
loc1 <- which(G_ID %in% First_NE)
Tumor <- R_data[TT,loc1]
Normal <- R_data[NN,loc1]

data_input <- as.matrix(rbind(Tumor,Normal))#data_input[1,]
colnames(data_input) <- as.character(First_NE)
label_input <- c(rep("tumor",length(which(Lable==1))),rep("normal",length(which(Lable==0))))#dim(label_input)
#subsets <- c(length(loc))#
subsets <- c(1:10,15,20,30,50)
ctrl <- rfeControl(functions=caretFuncs, 
                   method = "cv",repeats =2,number = 5,
                   returnResamp="final", verbose = TRUE)    ### 5-fold 交叉验证， 重复1 次
svm.ROC.Radial <- rfe(data_input, label_input, 
                      sizes=subsets,rfeControl=ctrl,method="svmRadial")
##  获取特征序后的结果
feature_ranks <- as.matrix(svm.ROC.Radial[["optVariables"]])
## 筛选前top 50 个特征
select <- feature_ranks[1:50]


## check 一下Base gene set 的分类效果
base_AUC <- t(as.matrix(unlist(Auc_all(Base_gene))))#0.981
## 把这些一起放进去看看分类效果
## 把选择的基因独自和base gene的效果加进去看看
AUC <- rbind()
for(ii in 1:length(select))
{
  all_gene <- union(Base_gene,select[ii])
   all_AUC <- t(as.matrix(unlist(Auc_all(all_gene))))
   AUC <- rbind(AUC,all_AUC[length(all_AUC)])
}
SS1 <- cbind(select[which(as.numeric(AUC)>base_AUC[length(base_AUC)])],
             AUC[which(as.numeric(AUC)>base_AUC[length(base_AUC)])])
## 对他们进行排序
SS11 <- SS1[order(as.numeric(SS1[,2]),decreasing = T),]
F_AUC <- rbind()
for(kk in 2:10)
{
Final <- union(Base_gene,SS11[1:kk,1])
Final_AUC1 <- t(as.matrix(unlist(Auc_all(Final))))
F_AUC <- rbind(F_AUC,c(kk,Final_AUC1[length(Final_AUC1)]))
}
## 筛选最大的即可
kkk1 <- which(as.numeric(F_AUC[,2]) == max(as.numeric(F_AUC[,2])))
#length(F_AUC[,2])
Module_gene <- Base_gene
for(jjj in 1:length(kkk1))
{
  Module_gene <- unique(union(Module_gene,SS11[1:(kkk1[jjj]+1),1]))
}

Module_search[[mmmm]] <- c(Module_gene,as.numeric(F_AUC[kkk1[1],2]))
AUCCC <- rbind(AUCCC,max(as.numeric(F_AUC[,2])))
###
#### 存储一下结果
#### 
#dir.create("Adding_12_27")
pat <- "./Adding_12_27/"
store_name <- paste(pat,as.character(DE_cur),sep="")
saveRDS(Module_search[[mmmm]], file = paste(store_name,"_module.RDS",sep=""))
}
#write.csv(AUC_store,paste(store_name,"_1_auc.csv"),quote=F,row.names=F)
### 最后再来module 去除重复即可
