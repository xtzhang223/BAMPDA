function [ drug_similarity,virus_similarity ] = integration_drug_virus_similarity(drug_gauss, virus_gauss)
%INT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

% wd = xlsread('weight_protein.xlsx');  %Ȩ�ؾ���
% wv = xlsread('weight_disease.xlsx');
% ss = xlsread('Protein_Sequence_Similarity_Matrix.xlsx'); %ҩ�ﻯѧ�ṹ�����Ծ���
% 
% vs = xlsread('Disease_Semantic_Similarity_Matrix.xlsx');          % ��������������
% 
% drug_similarity = wd.*ss+(~wd).*drug_gauss;     %ҩ�����������Ծ���
% 
% virus_similarity= wv.*vs+(~wv).*virus_gauss;  %�������������Ծ���




% wd = xlsread('weight_protein.xlsx');  %Ȩ�ؾ���
% wv = xlsread('weight_disease.xlsx');
% ss = xlsread('pp_decline.xlsx'); %ҩ�ﻯѧ�ṹ�����Ծ���
% 
% vs = xlsread('disease_decline.xlsx');          % ��������������
% 
% drug_similarity = wd.*ss+(~wd).*drug_gauss;     %ҩ�����������Ծ���
% 
% virus_similarity= wv.*vs+(~wv).*virus_gauss;  %�������������Ծ���




% wd = xlsread('weight_protein.xlsx');  %Ȩ�ؾ���
% wv = xlsread('weight_disease.xlsx');
ss = xlsread('Protein_Sequence_Similarity_Matrix.xlsx'); 

vs = xlsread('Disease_Semantic_Similarity_Matrix.xlsx');          

drug_similarity = ss+drug_gauss;     %ҩ�����������Ծ���

virus_similarity= vs+virus_gauss;  %�������������Ծ���
end

