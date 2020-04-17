clear; close all;clc;
dir = 'F:\bioinformatics\module_rank_20191209\20200319\Guide_of_mRank\';
Net = 'Regnetwork\';
Reg_dir = strcat(dir,Net);
Reg=importdata(strcat(Reg_dir,'renewed_RegNetwork_DMI_new.txt'));
%% 转换成 numeric
Reg_data =Reg.data;
Reg_data(:,3)=[];
S1 = Reg_data(:,1);
T1 = Reg_data(:,2);
W1= Reg_data(:,3);
%% 构建邻接矩阵
%    构建global 的转移矩阵
%     按一般的PageRank来做
node =unique(union(S1,T1));
Ad=zeros(length(node),length(node));
for(j =1:length(S1))
    l1 =find(S1(j)==node);
    l2 =find(T1(j)==node);
    if(~isempty(l1))
        if(~isempty(l2))
            Ad(l2,l1)=1;
        end
    end
end
%% 转移矩阵 正规化
%   列求和
No_Ad=Ad;
col_Sum = sum(Ad,1);
%length(find(col_Sum==0))
for(i =1:length(col_Sum))
    if(col_Sum>0)
        No_Ad(:,i)=No_Ad(:,i)/col_Sum(i);
    end
end
%% 读取module 信息
%Mod_dir ='module_integration\';
%Module = importdata(strcat(strcat(dir,Mod_dir),'reg_EX_module_all_095.txt'));
Module = importdata('module_final.csv');
% Module 转换成numeric
Module_used = Module.data;
Node_all = Module_used(:,:);
%Node_all(:,size(Node_all,2))=[];
PR_loc=zeros(size(Node_all,1),(size(Node_all,2)-1));
%% 构建local pagerank
for (i =1:size(Module_used,1))
      % i  = 1; 
                         k1 =  find(Module_used(i,:)>=1);
         %当前module的邻接矩阵
            loc_matrix = zeros(length(k1),length(k1));
            %当前module的node
        Node_inside = Module_used(i,k1);
        % 当前module的权重矩阵
          Loc_weight = zeros(length(k1),length(k1));
        % 在网络中查找所有的跟这些node相关的edge
                  for(j =1:length(Node_inside))
                      %j =3;
                          l1 = find(Node_inside(j)==S1);     
                           if(~isempty(l1))
                               W_in=W1(l1);
                               l2 = T1(l1);
                               N2 = setdiff(Node_inside,Node_inside(j));
                               for(kk = 1:length(N2))
                                   zz = find(N2(kk)==l2);
                                   if(~isempty(zz))
                                       loc1 = find(Node_inside(j)==Node_inside);
                                       loc2 = find(N2(kk) == Node_inside);
                                       loc_matrix(loc2,loc1)=1;
                                       Loc_weight(loc2,loc1)=W_in(zz);
                                   end
                               end
                          end
                  end
                 % 当前module的邻接矩阵正规化
                 % i---->j M(j,i)=1
                 % 就需要colsum
                 % 列求和
                 colSum = sum(loc_matrix);
                 nor_loc_matrix=loc_matrix;
                    for(nn =1:length(colSum))
                        if(colSum(nn)>0)
                        nor_loc_matrix(nn,:)=nor_loc_matrix(nn,:)/colSum;
                        end
                    end
                 % 有边就是有边，没有就是没有
                 % 当前module的邻接矩阵加权
                 W_loc_matrix=nor_loc_matrix.*Loc_weight;
                 % 当前module 内部的weighted PageRank算法
                               r = 0.85;
                 threshold = 1e-3;
                              N = length(W_loc_matrix);
                             PR = ( 1/N*ones(N,1));
                       restart = PR;
                            iter = 1;
                    delta_PR = Inf; 
                   while (delta_PR > threshold || iter>200)    
                            tic;
                            prev_PR = PR;               
                                     PR = r*W_loc_matrix* PR + (1-r)*restart;    
                            delta_PR = norm(PR-prev_PR);
                                  t(iter) = toc;
                                     iter = iter + 1;
                   end
                   %% 执行完local weighted PR
                   % 存储当前module的PR 值
                   for mm = 1:length(PR)
                   PR_loc(i,mm) = PR(mm);                       
                   end
end

% save('Modulepr.mat','PR_loc')
%save('Module_info.mat','Module_used')
% load('reg_results.mat')
% 创建超图
% 超图的邻接矩阵
Supergraph =zeros(size(Node_all,1),size(Node_all,1));
LS1=[];
count=1;
for(i =2115:size(Node_all,1))%当前module 和它下一个module

    for(j=1:size(Node_all,1))
        Z=0;
        if (i ~=j)
            oo1 = find(Node_all(i,:)>=1);
            oo2 =find(Node_all(j,:)>=1);
            M1_node =unique(Node_all(i,oo1));% module 1 
            M2_node= unique(Node_all(j,oo2));% module 2
             for(kk =1:length(M1_node)) % Module 1的所有node
                    kl1 = find(M1_node(kk)==node);
                    if(~isempty(kl1))%Module 2 的所有node
                      for(hh =1:length(M2_node))
                          kl2=find(M2_node(hh)==node);
                           if(~isempty(kl2))
                              Z= Z+PR_loc(i,kk)*Ad(kl2,kl1);%Ad： 网络的邻接矩阵
                           end
                      end
                    end
             end
               Supergraph(j,i)=Z;
               LS1(count,1)=i;
               LS1(count,2)=j;
               LS1(count,3)=Z;
               count=count+1;
        end
    end
end
%save()
length(find(LS1(:,3)>0))
% 3936 
mkdir('supergraph')
datacolumns = {'Source', 'Target', 'Weight'};
data = table( LS1(find(LS1(:,3)>0),1),...
                      LS1(find(LS1(:,3)>0),2),...
                      LS1(find(LS1(:,3)>0),3),...
                     'VariableNames', datacolumns); 
writetable(data ,strcat('supergraph\','reg_supergraph_relation.csv'));
%% 超图的PR
% 不加权重
L1 =sum(Supergraph,1);
LL=Supergraph;
% 需要对超图进行矩阵normalize 不然运行PageRank不会收敛
      for(j=1:length(L1))
          if(L1(j)>0)
          LL(:,j)=LL(:,j)/(L1(j));
          end
      end
                    r = 0.85;
      threshold = 1e-5;
                  N = length(LL);
                PR = (1/N*ones(N,1));
         restart = PR;
              iter = 1;
     delta_PR = Inf; 
         while (delta_PR > threshold || iter<200)    
                            tic;
                         prev_PR = PR;               
                                 PR = r*LL* PR + (1-r)*restart;    
                        delta_PR= norm(PR-prev_PR);
                             t(iter)=toc;
                               iter = iter + 1;
         end
        PR_global=PR;
        % 由于筛选的时候，以Auc 0.97为阈值
        % 故此次排序的时候只考虑Global PR
        Final_results = [Module_used,PR_global];
        size(Final_results)
datacolumns ={'V1','V2','V3','V4','V5',...
                         'V6','V7','V8','V9','V10',...
                         'V11','V12','V13','V14','V15',...
                         'V16','V17','V18'};
data_final = table( Final_results(:,1),Final_results(:,2),Final_results(:,3),Final_results(:,4),Final_results(:,5),...
                              Final_results(:,6),Final_results(:,7),Final_results(:,8),Final_results(:,9),Final_results(:,10),...
                              Final_results(:,11),Final_results(:,12),Final_results(:,13),Final_results(:,14),Final_results(:,15), ...
                              Final_results(:,16),Final_results(:,17),Final_results(:,18), ...
                     'VariableNames', datacolumns); 
writetable( data_final ,strcat('supergraph\','reg_supergraph_info_results.csv'));


% 对各个指标单独排序
% save('Final_results.mat','Final_results','-v6')
P_Rank = sort(unique(PR_global),'descend');
A1 = [];
for(i = 1:size(Module_used,1))
%     i =1
    k1 = find(Module_used(i,:)<1&Module_used(i,:)>0);
    A1(i,1)=Module_used(i,k1);
end
A_Rank = sort(unique(A1),'descend');
LABEL = 1:size(Module_used,1);
F_store=[LABEL',A1,PR_global];
%D_Rank = sort(unique(Final_results(:,2)),'descend');
% 整合各个指标的排序
% 
Rank=[];
for(i =1:size(F_store,1))
   
    k1 = find(A1(i) == A_Rank);
    k2 = find(PR_global(i) == P_Rank);
    %k3 = find(Final_results(i,2) == D_Rank);
    kk=[k1,k2];
    k4 = sum(kk)/2;
    Rank(i,1)=k1;
    Rank(i,2)=k2;
    %Rank(i,3)=k3;
    Rank(i,3)=k4;
end
% size(data_final)
 datacolumns1 = { 'Lable',...
                               'Auc','A_Rank','PR','P_Rank',...
                              'F_rank'};
 data_f= table(F_store(:,1),F_store(:,2),Rank(:,1),...
                        F_store(:,3),Rank(:,2),...
                        Rank(:,3),...
                     'VariableNames', datacolumns1); 
%                  size(data_f)
     
writetable(data_f ,strcat('supergraph\','reg_supergraph_final_results_with_rank.csv'));
% save reg_results.mat
% 超图中每个super graph 的size
Store_size =zeros(size(Node_all,1),2);
% size(Node_all,2)
for(i = 1:size(Node_all,1))
    k1 = find(Node_all(i,:)>=1);
    nodes =length(k1);
    Store_size(i,1)=i;
    Store_size(i,2)=nodes;
end
datacolumns1 = {'Super_graph',  'SIZE'};
data_f= table(Store_size(:,1),Store_size(:,2), 'VariableNames', datacolumns1);
writetable(data_f ,strcat('supergraph\','reg_supergraph_size.csv'));

%%最后整合成一个大表

datacolumns ={'V1','V2','V3','V4','V5',...
                         'V6','V7','V8','V9','V10',...
                         'V11','V12','V13','V14','V15',...
                         'V16','V17','V18',...
                         'Auc','A_Rank','PR','P_Rank','F_rank', 'SIZE'};
data_final2 = table( Final_results(:,1),Final_results(:,2),Final_results(:,3),Final_results(:,4),Final_results(:,5),...
                              Final_results(:,6),Final_results(:,7),Final_results(:,8),Final_results(:,9),Final_results(:,10),...
                              Final_results(:,11),Final_results(:,12),Final_results(:,13),Final_results(:,14),Final_results(:,15), ...
                              Final_results(:,16),Final_results(:,17),Final_results(:,18), ...
                              F_store(:,2),Rank(:,1),F_store(:,3),Rank(:,2),Rank(:,3),...
                             Store_size(:,2),...
                            'VariableNames', datacolumns); 
writetable(data_final2 ,strcat('supergraph\','reg_supergraph_store.csv'));

save('Block_rank.mat')