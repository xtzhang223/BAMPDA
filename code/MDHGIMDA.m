function [outputArg1] = MDHGIMDA(interaction_matrix1)
%MDHGIMDA 此处显示有关此函数的摘要
%   此处显示详细说明
alpha1=0.1;                                                              % define the parameters
% alpha1=0.5;
maxit=10e-06;
alpha=0.4;
Adaj_dm=interaction_matrix1;
AAdaj_dm=Adaj_dm;
[X,E] = lrr(AAdaj_dm,alpha1);
AAdaj_dm= AAdaj_dm*X;
change=1;
ADAJ_M= AAdaj_dm';
ADAJ_DM= AAdaj_dm;
ADAJ_sco=AAdaj_dm;
drug_gauss = similarity_drug(AAdaj_dm);
microbe_gauss = similarity_microbe(AAdaj_dm);
[drug_integration_similarity,microbe_integration_similarity] = integration_drug_virus_similarity(drug_gauss, microbe_gauss ); %调用函数获得整合药物副作用，药物化学结构相似性，高斯相似性的药物相似性矩阵，病毒相似性也是如此
ID=drug_integration_similarity;
IM=microbe_integration_similarity;
%计算得分
change=1;
Norm_adaj_dm1=ID;
Norm_adaj_m1=IM;
d1=sum( Norm_adaj_dm1,2);
d2=sum( Norm_adaj_m1,2);
D1=diag(1./sqrt(d1));
D1(D1==Inf)=0;
D2=diag(1./sqrt(d2));
D2(D2==Inf)=0;
Norm_adaj_dm2=D1* Norm_adaj_dm1*D1;                              % normalization the integrated similarity of KD
Norm_adaj_m2=D2* Norm_adaj_m1*D2;                                % normalization the integrated similarity of KM
while( change>maxit)
    NEW_adaj_dm=alpha* Norm_adaj_dm2* ADAJ_DM*Norm_adaj_m2+(1-alpha)* AAdaj_dm;          %the matrix of predicting score between disease-miRNA
    %      NEW_adaj_m=alpha* Norm_adaj_m2* ADAJ_M+(1-alpha)* AAdaj_dm';            %the matrix of predicting score between miRNA-lncRNA
    New_adaj_sco=NEW_adaj_dm;                         %predicted score of miRNA-disease associations
    change=sum(sum((abs(New_adaj_sco- ADAJ_sco)))) ;                      % the iterative change
    ADAJ_DM=NEW_adaj_dm;                                                  %assign NEW_adaj_dm to ADAJ_DM
    %      ADAJ_M= NEW_adaj_m;                                                   %assign NEW_adaj_m to ADAJ_M
    ADAJ_sco=ADAJ_DM;
end
outputArg1 = ADAJ_sco;
end

