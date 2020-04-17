rm(list=ls())
setwd("F:/bioinformatics/module_rank_20191209/Results/m_Rank_Top_10_validated_with_other_datasest")
ID_used <- as.numeric(as.matrix(read.csv("HCC_biomarkers_detected.csv",header = T))[,1])
#dim(ID_used)#ID_used[1,]
##读取下载之后的数据
traindata <- as.matrix(read.csv("TCGA_lihc_origin_unique.csv",header = T))
#dim(traindata)#traindata[1,1]
used_T <- traindata[,-1]
rownames(used_T)<- traindata[,1]
library(org.Hs.eg.db)
T1 <- rownames(used_T)
T_ID <- mapIds(org.Hs.eg.db, keys=T1, 
               column="ENTREZID",keytype="SYMBOL",
               multiVals="first")
kss <- which(T_ID != "")#length(kss)
used_T <- as.matrix(used_T[kss,])#dim(used_T)
trainlabel <- as.matrix(read.csv("gene_data_lable.csv",header=T)) #训练数据 sample 标签#dim(trainlabel)
trainlabel_new <- trainlabel#trainlabel[1]
#T_ID <- as.numeric(traindata[,1])# 训练数据gene ID#length(which(D_ID %in% T_ID))#length(which(T_ID %in% D_ID))
T_ID <- as.numeric(T_ID[kss]) 
rownames(used_T)<- T_ID
length(T_ID)
### 读取DE 信息
DE  <- as.matrix(read.table("DE_gene_with_ID.txt",sep="\t",header = F))[,3]
###45436
Da_22 <- as.matrix(read.csv("arraydata22058.csv",header = T))
Da_used_22 <-as.matrix(Da_22[,-1])#dim(Da_used_22)
D_ID_22 <-as.numeric(Da_22[,1])#length(D_ID_22)#13560
rownames(Da_used_22) <- D_ID_22
#rownames(Da_used_22)[1:2]
Inter <- intersect(as.numeric(rownames(Da_used_22)),
                   as.numeric(rownames(used_T)))#length(Inter)
##
used_T_In <- used_T[which(T_ID %in% Inter),] 
#dim(used_T_In)#rownames(used_T_In)[1]
#T_ID[which(T_ID %in% Inter)[1]]
used_22 <- Da_used_22[which(D_ID_22 %in% Inter),]
lable_22 <- apply(as.matrix(colnames(Da_22)[-1]),1, function(x) substring(x,2,2))
lable_22_new <-lable_22
lable_22_new[which(as.numeric(lable_22)==0)] <-c("nontumor")
lable_22_new[which(as.numeric(lable_22)==1)] <-c("tumor")
lable_22_new <- as.factor(lable_22_new)
used_Module <- ID_used
used_Module_DE <- intersect(used_Module,as.numeric(DE))
gg <- which(as.numeric(rownames(used_22)) %in% used_Module_DE)
length(gg)
library(pracma)
k1 <- which(as.numeric(rownames(used_22)) %in% used_Module)
k2 <- which(as.numeric(rownames(used_T_In)) 
            %in% as.numeric(rownames(used_22))[k1])#k2
#length(k2)
#T_data <- as.matrix(traindata[,3:dim(traindata)[2]])
T_data <- as.matrix(used_T_In)
T_data <- as.matrix(apply(T_data,2,function(x)as.numeric(x)))
z11 <- max(T_data)
z12 <- min(T_data)
T_data <- apply(T_data,2,function(x)(x-z12)/(z11-z12))
#max(T_data)#min(T_data)
T_used <- as.matrix(T_data[k2,])
trainlabel_new[which(as.numeric(trainlabel)==1)] <-c("tumor")
trainlabel_new[which(as.numeric(trainlabel)==0)] <-c("Nontumor")
trainlabel_new <- as.factor(trainlabel_new)
train_used <- data.frame(t(T_used),trainlabel_new)#dim(train_used)
x <-t(T_used)# dim(x)
y <-trainlabel_new
set.seed(1234)
library("e1071")
library(pROC)
svm_tune <- tune.svm(x,y,kernel="radial",
                     cost=10^(-1:4), gamma=c(.5,1,2))#print(svm_tune)#mode(svm_tune)
best_mod <- svm_tune$best.model
pred_after <- predict(best_mod,x,decision.values=TRUE)
f<-attr(pred_after,"decision.values")#提取决策值
prb1<-1/(1+exp(-f))#把决策值使用sigmoid 函数映射到0~1之间
show1 <- roc(as.vector(trainlabel),as.vector(prb1),plot=F,print.auc=F,print.thres=F,col="red",lty=1)
auc_train <-as.numeric(show1[["auc"]])
auc_train
Data_input_22 <- as.matrix(apply(as.matrix(used_22),2,
                                 function(x)as.numeric(x)))
#a11 <-max(Data_input_22) 
#a12 <- min(Data_input_22)
Data_input_22 <- apply(Data_input_22,2,function(x)(x-z12)/(z11-z12))
dim(Data_input_22)
Data_input_22 <- t(Data_input_22[k1,])
dim(Data_input_22)
#Data_input_22 <- apply(Data_input_22,1,function(x) as.numeric(x))#dim(Data_input)
pre_new_22 <- predict(best_mod,Data_input_22,decision.values=TRUE) 
f_new_22 <- attr(pre_new_22,"decision.values")#提取决策值
prb1_new_22 <- 1/(1+exp(-f_new_22))#

show_22<-roc(as.vector(lable_22_new),as.vector(prb1_new_22),
             plot=F,print.auc=F,print.thres=F,col="blue",
             lty=1,add=F) #
auc_22 <-as.numeric(show_22[["auc"]])
auc_22
zz <- show_22[["sensitivities"]]+show_22[["specificities"]]
k <- which(zz==max(zz))
se_cur <- show_22[["sensitivities"]][k[1]]
sp_cur <- show_22[["specificities"]][k[1]]
auc_cur <- as.numeric(show_22[["auc"]])
acc_cur <- (sp_cur+se_cur)/2
pre_cur <- se_cur/(se_cur+1-sp_cur)
f1_cur <- 2*pre_cur*se_cur/(pre_cur+se_cur)
Index <- c(se_cur,sp_cur,acc_cur,f1_cur,auc_cur)
Index
saveRDS(show_22,"GSE22058_0406_new.RDS")
