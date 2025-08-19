A=load('Protein_Disease_adj.txt');                                                       %注意  RBM 代表 MDHGI, INMB 代表 BNNR
nd=max(A(:,1));                                                           % nd:the number of drug
nv=max(A(:,2));                                                           % nv:the number of virus
[pp,qq]=size(A);                                                          % pp:the number of known drug-virus associations,pp=933,qq=3
% adaj_dm=zeros(nd,nv);
% for i=1:pp
%         adaj_dm(A(i,1),A(i,2))=1;                                           %adaj_dm: adajency matrix for the drug-virus association network,已知联系的置1
% end
interaction_matrix = xlsread('Protein_Disease_Associations.xlsx');
[row,col]=size(interaction_matrix);
number=100;                                                               %define the parameters
C1_BRM=[];
C1_INBM=[];
C1_ensemble1=[];
C1_ensemble2=[];
C1_ensemble3=[];
C1_ensemble4=[];
C1_ensemble5=[];
C1_ensemble6=[];
C1_ensemble7=[];
C1_ensemble8=[];
C1_ensemble9=[];
five_k_RBM=zeros(number,pp);
five_k_INBM=zeros(number,pp);
five_k_ensemble1=zeros(number,pp);
five_k_ensemble2=zeros(number,pp);
five_k_ensemble3=zeros(number,pp);
five_k_ensemble4=zeros(number,pp);
five_k_ensemble5=zeros(number,pp);
five_k_ensemble6=zeros(number,pp);
five_k_ensemble7=zeros(number,pp);
five_k_ensemble8=zeros(number,pp);
five_k_ensemble9=zeros(number,pp);
% for n=1:9
% C1_ensemble=zeros(9,pp);
% end
% for n=1:9
% C2=zeros(9,pp);
% end

% BNNR 模型的一些参数
maxiter = 300;
alpha = 1;
beta = 10;
tol1 = 2*1e-3;
tol2 = 1*1e-5;


array1= xlsread('C:\Users\20379\Desktop\paint-roc\array_P-D.xlsx');
tic;
for k=58:number
    fprintf('第%d轮\n',k);
    current_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
    disp(current_time);
    array=array1(k,:);                                      %array是100*933的矩阵，每一行是1-933的随机序列
    for circle=0:4
                                                  %因为933除以5除不尽，所以前四折是187，最后一折为185
            nu=181;
            interaction_matrix1=interaction_matrix;
            interaction_matrix1_index=logical(zeros(row,col)~=0);    %interaction_matrix1_index为训练样本对应位置的矩阵
            for i=1:pp
                interaction_matrix1_index(A(i,1),A(i,2))=1;
            end
            new_array = array(1,1+circle*181:(circle+1)*181);      % 933个药物-病毒关联对，前面四个nu=187，最后一个n=185
            
            for j=1:181                                            %j=1:nu
                o=A(new_array(1,j),1);
                l=A(new_array(1,j),2);
                interaction_matrix1(o,l)=0;                       % 依次置将1变为0
                interaction_matrix1_index(o,l)=0;
            end
            
            drug_gauss = xlsread('Protein_Gaussian_Similarity_Matrix.xlsx');
            virus_gauss = xlsread('Disease_Gaussian_Similarity_Matrix.xlsx');
            [drug_integration_similarity,virus_integration_similarity] = integration_drug_virus_similarity(drug_gauss, virus_gauss ); %调用函数获得整合药物副作用，药物化学结构相似性，高斯相似性的药物相似性矩阵，病毒相似性也是如此
            
            %调用三种方法获取得分矩阵
%             [predict_score_matrix1,numhid]=RBM_model(interaction_matrix1);   %调用RBM方法，numhid是number of hidden layers
%             [predict_score_matrix_drug,threshold]=integrated_neighbor_model_drug( interaction_matrix1,drug_integration_similarity);   %整合的领域的方法,用药物-病毒相似矩阵和药物相似矩阵
%             [predict_score_matrix_virus,threshold]=integrated_neighbor_model_virus( interaction_matrix1,virus_integration_similarity);
%             predict_score_matrix2=(predict_score_matrix_drug+predict_score_matrix_virus)/2;
            
            % MDHGIMDA 预测的关联得分
            predict_score_matrix1 = MDHGIMDA(interaction_matrix1);
            
            %BNNR 预测的关联得分
            Wdd = drug_integration_similarity;
            Wvv= virus_integration_similarity;
            Wvd = interaction_matrix1';
            [dn,dr] = size(Wvd);

            T = [Wdd, Wvd'; Wvd, Wvv];
            [t1, t2] = size(T);
            trIndex = double(T ~= 0);
            [WW,iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
            M_recovery = WW((t1-dn+1) : t1, 1 : dr);
            predict_score_matrix2 = M_recovery';
            
            predict_score_matrix3_1=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.1);     %调用集成的方法
            predict_score_matrix3_2=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.2);     %调用集成的方法
            predict_score_matrix3_3=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.3);     %调用集成的方法
            predict_score_matrix3_4=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.4);     %调用集成的方法
            predict_score_matrix3_5=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.5);     %调用集成的方法
            predict_score_matrix3_6=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.6);     %调用集成的方法
            predict_score_matrix3_7=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.7);     %调用集成的方法
            predict_score_matrix3_8=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.8);     %调用集成的方法
            predict_score_matrix3_9=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.9);     %调用集成的方法
            %     predict_score_matrix3=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index);     %调用集成的方法
            %
            %第一种方法
            Sco1= predict_score_matrix1;
            final_score=Sco1( interaction_matrix==0);
            for i=1:181
                q=A(new_array(1,i),1);
                w=A(new_array(1,i),2);
                s_score=Sco1(q,w);
                T=[s_score;final_score];                                                %sort the score
                index1=sort(T,'descend');
                index1_1=find(index1==s_score);                                          %find the sort of changed variate
                if length(index1_1)~=1
                    index1_2=mean(index1_1);                                              %decide the number of the same value
                else
                    index1_2=index1_1;                                                    %get the average rank between same value
                end
                five_g1(i)=index1_2;
            end
            
            %第二种：INBM方法
            Sco2= predict_score_matrix2;
            final_score=Sco2( interaction_matrix==0);
            for i=1:181
                q=A(new_array(1,i),1);
                w=A(new_array(1,i),2);
                s_score=Sco2(q,w);
                T=[s_score;final_score];                                                %sort the score
                index2=sort(T,'descend');
                index2_1=find(index2==s_score);                                          %find the sort of changed variate
                if length(index2_1)~=1
                    index2_2=mean(index2_1);                                              %decide the number of the same value
                else
                    index2_2=index2_1;                                                    %get the average rank between same value
                end
                five_g2(i)=index2_2;
            end
            
            %第三种方法：ensemble
            
            %        Sco3= predict_score_matrix3;
            %        final_score=Sco3( interaction_matrix==0);
            %        for i=1:187
            %          q=A(new_array(1,i),1);
            %          w=A(new_array(1,i),2);
            %          s_score=Sco3(q,w);
            %          T=[s_score;final_score];                                                %sort the score
            %          index3=sort(T,'descend');
            %          index3_1=find(index3==s_score);                                          %find the sort of changed variate
            %          if length(index3_1)~=1
            %            index3_2=mean(index3_1);                                              %decide the number of the same value
            %          else
            %            index3_2=index3_1;                                                    %get the average rank between same value
            %         end
            %         five_g3(i)=index3_2;
            %        end
            
            five1 = conputing_ensemble_ranking( predict_score_matrix3_1,interaction_matrix,181,A,new_array); %conputing_ensemble_ranking 计算集成算法的排名
            five2 = conputing_ensemble_ranking( predict_score_matrix3_2,interaction_matrix,181,A,new_array);
            five3 = conputing_ensemble_ranking( predict_score_matrix3_3,interaction_matrix,181,A,new_array);
            five4 = conputing_ensemble_ranking( predict_score_matrix3_4,interaction_matrix,181,A,new_array);
            five5 = conputing_ensemble_ranking( predict_score_matrix3_5,interaction_matrix,181,A,new_array);
            five6 = conputing_ensemble_ranking( predict_score_matrix3_6,interaction_matrix,181,A,new_array);
            five7 = conputing_ensemble_ranking( predict_score_matrix3_7,interaction_matrix,181,A,new_array);
            five8 = conputing_ensemble_ranking( predict_score_matrix3_8,interaction_matrix,181,A,new_array);
            five9 = conputing_ensemble_ranking( predict_score_matrix3_9,interaction_matrix,181,A,new_array);
            
            C1_BRM=[C1_BRM,five_g1];
            C1_INBM=[C1_INBM,five_g2];
            
            C1_ensemble1=[C1_ensemble1,five1];
            C1_ensemble2=[C1_ensemble2,five2];
            C1_ensemble3=[C1_ensemble3,five3];
            C1_ensemble4=[C1_ensemble4,five4];
            C1_ensemble5=[C1_ensemble5,five5];
            C1_ensemble6=[C1_ensemble6,five6];
            C1_ensemble7=[C1_ensemble7,five7];
            C1_ensemble8=[C1_ensemble8,five8];
            C1_ensemble9=[C1_ensemble9,five9];
            
        end
        
        if circle==4
            nu=181;
            interaction_matrix1=interaction_matrix;
            interaction_matrix1_index=logical(zeros(row,col)~=0);
            for i=1:pp
                interaction_matrix1_index(A(i,1),A(i,2))=1;
            end
            new_array = array(1,1+circle*181:pp);
            for j=1:181
                o=A(new_array(1,j),1);
                l=A(new_array(1,j),2);
                interaction_matrix1(o,l)=0;
                interaction_matrix1_index(o,l)=0;
            end
            drug_gauss = similarity_drug(interaction_matrix1);
            virus_gauss = similarity_microbe(interaction_matrix1);
            [drug_integration_similarity,virus_integration_similarity] = integration_drug_virus_similarity(drug_gauss, virus_gauss ); %调用函数获得整合药物副作用，药物化学结构相似性，高斯相似性的药物相似性矩阵，病毒相似性也是如此
            
            %调用三种方法获取得分矩阵
%             [predict_score_matrix1,numhid]=RBM_model(interaction_matrix1);   %调用RBM方法，numhid是number of hidden layers
%             [predict_score_matrix_drug,threshold]=integrated_neighbor_model_drug( interaction_matrix1,drug_integration_similarity);   %整合的领域的方法,用药物-病毒相似矩阵和药物相似矩阵
%             [predict_score_matrix_virus,threshold]=integrated_neighbor_model_virus( interaction_matrix1,virus_integration_similarity);
%             predict_score_matrix2=(predict_score_matrix_drug+predict_score_matrix_virus)/2;
            
            predict_score_matrix1 = MDHGIMDA(interaction_matrix1);
            
            %BNNR 预测的关联得分
            Wdd = drug_integration_similarity;
            Wvv= virus_integration_similarity;
            Wvd = interaction_matrix1';
            [dn,dr] = size(Wvd);

            T = [Wdd, Wvd'; Wvd, Wvv];
            [t1, t2] = size(T);
            trIndex = double(T ~= 0);
            [WW,iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
            M_recovery = WW((t1-dn+1) : t1, 1 : dr);
            predict_score_matrix2 = M_recovery';
            
            predict_score_matrix3_1=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.1);     %调用集成的方法
            predict_score_matrix3_2=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.2);     %调用集成的方法
            predict_score_matrix3_3=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.3);     %调用集成的方法
            predict_score_matrix3_4=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.4);     %调用集成的方法
            predict_score_matrix3_5=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.5);     %调用集成的方法
            predict_score_matrix3_6=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.6);     %调用集成的方法
            predict_score_matrix3_7=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.7);     %调用集成的方法
            predict_score_matrix3_8=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.8);     %调用集成的方法
            predict_score_matrix3_9=ensemble_average(predict_score_matrix1,predict_score_matrix2,interaction_matrix1_index,0.9);
            
            
            
            %调用第一种方法RBM
            Sco1= predict_score_matrix1;
            final_score=Sco1( interaction_matrix==0);
            for i=1:181
                q=A(new_array(1,i),1);
                w=A(new_array(1,i),2);
                s_score=Sco1(q,w);
                T=[s_score;final_score];                                                %sort the score
                index1=sort(T,'descend');
                index1_1=find(index1==s_score);                                          %find the sort of changed variate
                
                if length(index1_1)~=1
                    index1_2=mean(index1_1);                                              %decide the number of the same value
                else
                    index1_2=index1_1;                                                    %get the average rank between same value
                end
                five_g1_1(i)=index1_2;
            end
            
            
            %第二种：INBM方法
            Sco2= predict_score_matrix2;
            final_score=Sco2( interaction_matrix==0);
            for i=1:181
                q=A(new_array(1,i),1);
                w=A(new_array(1,i),2);
                s_score=Sco2(q,w);
                T=[s_score;final_score];                                                %sort the score
                index2=sort(T,'descend');
                index2_1=find(index2==s_score);                                          %find the sort of changed variate
                if length(index2_1)~=1
                    index2_2=mean(index2_1);                                              %decide the number of the same value
                else
                    index2_2=index2_1;                                                    %get the average rank between same value
                end
                five_g2_1(i)=index2_2;
            end
            
            
            %         %第三种方法：ensemble
            %         Sco3= predict_score_matrix3(n);
            %         final_score=Sco3( interaction_matrix==0);
            %         for i=1:185
            %           q=A(new_array(1,i),1);
            %           w=A(new_array(1,i),2);
            %           s_score=Sco3(q,w);
            %           T=[s_score;final_score];                                                %sort the score
            %           index3=sort(T,'descend');
            %           index3_1=find(index3==s_score);                                          %find the sort of changed variate
            %           if length(index3_1)~=1
            %             index3_2=mean(index3_1);                                              %decide the number of the same value
            %           else
            %             index3_2=index3_1;                                                    %get the average rank between same value
            %           end
            %           five_g3_1(i)=index3_2;
            %         end
            
            
            five1 = conputing_ensemble_ranking( predict_score_matrix3_1,interaction_matrix,181,A,new_array); %conputing_ensemble_ranking 计算集成算法的排名
            five2 = conputing_ensemble_ranking( predict_score_matrix3_2,interaction_matrix,181,A,new_array);
            five3 = conputing_ensemble_ranking( predict_score_matrix3_3,interaction_matrix,181,A,new_array);
            five4 = conputing_ensemble_ranking( predict_score_matrix3_4,interaction_matrix,181,A,new_array);
            five5 = conputing_ensemble_ranking( predict_score_matrix3_5,interaction_matrix,181,A,new_array);
            five6 = conputing_ensemble_ranking( predict_score_matrix3_6,interaction_matrix,181,A,new_array);
            five7 = conputing_ensemble_ranking( predict_score_matrix3_7,interaction_matrix,181,A,new_array);
            five8 = conputing_ensemble_ranking( predict_score_matrix3_8,interaction_matrix,181,A,new_array);
            five9 = conputing_ensemble_ranking( predict_score_matrix3_9,interaction_matrix,181,A,new_array);
            
            C1_BRM=[C1_BRM,five_g1_1];
            C1_INBM=[C1_INBM,five_g2_1];
            
            C1_ensemble1=[C1_ensemble1,five1];
            C1_ensemble2=[C1_ensemble2,five2];
            C1_ensemble3=[C1_ensemble3,five3];
            C1_ensemble4=[C1_ensemble4,five4];
            C1_ensemble5=[C1_ensemble5,five5];
            C1_ensemble6=[C1_ensemble6,five6];
            C1_ensemble7=[C1_ensemble7,five7];
            C1_ensemble8=[C1_ensemble8,five8];
            C1_ensemble9=[C1_ensemble9,five9];
        end
    end
    
    
    
    
    
    
    
    
    
    five_k_RBM(k,:)=C1_BRM;
    five_k_INBM(k,:)=C1_INBM;
    C1_BRM=[];
    C1_INBM=[];
    five_k_ensemble1(k,:)=C1_ensemble1;
    five_k_ensemble2(k,:)=C1_ensemble2;
    five_k_ensemble3(k,:)=C1_ensemble3;
    five_k_ensemble4(k,:)=C1_ensemble4;
    five_k_ensemble5(k,:)=C1_ensemble5;
    five_k_ensemble6(k,:)=C1_ensemble6;
    five_k_ensemble7(k,:)=C1_ensemble7;
    five_k_ensemble8(k,:)=C1_ensemble8;
    five_k_ensemble9(k,:)=C1_ensemble9;
    C1_ensemble1=[];
    C1_ensemble2=[];
    C1_ensemble3=[];
    C1_ensemble4=[];
    C1_ensemble5=[];
    C1_ensemble6=[];
    C1_ensemble7=[];
    C1_ensemble8=[];
    C1_ensemble9=[];
    
    
    % five_k_ensemble_1=C2(1,:);
    % five_k_ensemble_2=C2(2,:);
    % five_k_ensemble_3=C2(3,:);
    % five_k_ensemble_4=C2(4,:);
    % five_k_ensemble_5=C2(5,:);
    % five_k_ensemble_6=C2(6,:);
    % five_k_ensemble_7=C2(7,:);
    % five_k_ensemble_8=C2(8,:);
    % five_k_ensemble_9=C2(9,:);
    
    %     for n=1:9
    %     C2(i,:)=[];
    %     end
    %
    
end
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\result_58-100\five_fold_MDHGI.xlsx',five_k_RBM);
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\result_58-100\five_fold_BNNR.xlsx',five_k_INBM);

xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\result_58-100\five_fold_ensemble1.xlsx',five_k_ensemble1);
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\result_58-100\five_fold_ensemble2.xlsx',five_k_ensemble2);
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\result_58-100\five_fold_ensemble3.xlsx',five_k_ensemble3);
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\result_58-100\five_fold_ensemble4.xlsx',five_k_ensemble4);
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\result_58-100\five_fold_ensemble5.xlsx',five_k_ensemble5);
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\result_58-100\five_fold_ensemble6.xlsx',five_k_ensemble6);
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\result_58-100\five_fold_ensemble7.xlsx',five_k_ensemble7);
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\result_58-100\five_fold_ensemble8.xlsx',five_k_ensemble8);
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\result_58-100\five_fold_ensemble9.xlsx',five_k_ensemble9);


toc;
fprintf('five_fold结束\n');


