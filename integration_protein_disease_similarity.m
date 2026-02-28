function [ protein_similarity,disease_similarity ] = integration_protein_disease_similarity(protein_gauss, disease_gauss)
%INT 此处显示有关此函数的摘要
%   此处显示详细说明











% wd = xlsread('weight_protein.xlsx');  %权重矩阵
% wv = xlsrea'weight_disease.xlsx');
ss = xlsread('Protein_Sequence_Similarity_Matrix.xlsx'); 

vs = xlsread('Disease_Semantic_Similarity_Matrix.xlsx');          

protein_similarity = ss+protein_gauss;     %药物最终相似性矩阵

disease_similarity= vs+disease_gauss;  %病毒最终相似性矩阵
end

