
position_BNNR=xlsread('five_fold_ensemble3.xlsx');

%已知关联矩阵
interaction=xlsread('Protein_Disease_Associations.xlsx');

%已知关联对
sID=xlsread('Protein_Disease_adj.xlsx');

overallauc_BNNR=auc(position_BNNR);


auc_mean_BNNR=mean(overallauc_BNNR)

