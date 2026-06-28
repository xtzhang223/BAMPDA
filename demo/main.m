clear;
clc;

%% 路径设置
result_dir = 'C:\Users\zxt\Desktop\run_MDHGIBNNR\demo';


%% ============================================================

A = load(fullfile(result_dir, 'Protein_Disease_adj.txt'));                
nd = max(A(:,1));
nv = max(A(:,2));
[pp,qq] = size(A);

interaction_matrix = xlsread(fullfile(result_dir, 'Protein_Disease_Associations.xlsx'));
[row,col] = size(interaction_matrix);
number = 100;                                                             % define the parameters

C1_ensemble3 = [];
five_k_ensemble3 = zeros(number,pp);

% BNNR 模型的一些参数
maxiter = 300;
alpha = 1;
beta = 10;
tol1 = 2*1e-3;
tol2 = 1*1e-5;

array1 = xlsread(fullfile(result_dir, 'array_P-D.xlsx'));

tic;
for k = 95:number
    fprintf('第%d轮\n',k);
    current_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
    disp(current_time);
    array = array1(k,:);                                  

    for circle = 0:4
        if circle < 4
            nu = 181;
            interaction_matrix1 = interaction_matrix;
            interaction_matrix1_index = logical(zeros(row,col) ~= 0);
            for i = 1:pp
                interaction_matrix1_index(A(i,1),A(i,2)) = 1;
            end
            new_array = array(1,1+circle*181:(circle+1)*181);

            for j = 1:181
                o = A(new_array(1,j),1);
                l = A(new_array(1,j),2);
                interaction_matrix1(o,l) = 0;
                interaction_matrix1_index(o,l) = 0;
            end

            protein_gauss = xlsread(fullfile(result_dir, 'Protein_Gaussian_Similarity_Matrix.xlsx'));
            disease_gauss = xlsread(fullfile(result_dir, 'Disease_Gaussian_Similarity_Matrix.xlsx'));
            [protein_integration_similarity,disease_integration_similarity] = ...
                integration_protein_disease_similarity(protein_gauss, disease_gauss);

            % MDHGIMDA 预测的关联得分
            predict_score_matrix1 = MDHGIMDA(interaction_matrix1);

            % BNNR 预测的关联得分
            Wdd = protein_integration_similarity;
            Wvv = disease_integration_similarity;
            Wvd = interaction_matrix1';
            [dn,dr] = size(Wvd);

            T = [Wdd, Wvd'; Wvd, Wvv];
            [t1,t2] = size(T);
            trIndex = double(T ~= 0);
            [WW,iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
            M_recovery = WW((t1-dn+1) : t1, 1 : dr);
            predict_score_matrix2 = M_recovery';

            predict_score_matrix3_3 = ensemble_average( ...
                predict_score_matrix1, ...
                predict_score_matrix2, ...
                interaction_matrix1_index, ...
                0.3);

            % 第一种方法
            Sco1 = predict_score_matrix1;
            final_score = Sco1(interaction_matrix == 0);
            for i = 1:181
                q = A(new_array(1,i),1);
                w = A(new_array(1,i),2);
                s_score = Sco1(q,w);
                T = [s_score;final_score];
                index1 = sort(T,'descend');
                index1_1 = find(index1 == s_score);
                if length(index1_1) ~= 1
                    index1_2 = mean(index1_1);
                else
                    index1_2 = index1_1;
                end
                five_g1(i) = index1_2;
            end

            % 第二种：INBM 方法
            Sco2 = predict_score_matrix2;
            final_score = Sco2(interaction_matrix == 0);
            for i = 1:181
                q = A(new_array(1,i),1);
                w = A(new_array(1,i),2);
                s_score = Sco2(q,w);
                T = [s_score;final_score];
                index2 = sort(T,'descend');
                index2_1 = find(index2 == s_score);
                if length(index2_1) ~= 1
                    index2_2 = mean(index2_1);
                else
                    index2_2 = index2_1;
                end
                five_g2(i) = index2_2;
            end

            five3 = conputing_ensemble_ranking( ...
                predict_score_matrix3_3, ...
                interaction_matrix, ...
                181, ...
                A, ...
                new_array);

            C1_ensemble3 = [C1_ensemble3,five3];
        end

        if circle == 4
            nu = 181;
            interaction_matrix1 = interaction_matrix;
            interaction_matrix1_index = logical(zeros(row,col) ~= 0);
            for i = 1:pp
                interaction_matrix1_index(A(i,1),A(i,2)) = 1;
            end
            new_array = array(1,1+circle*181:pp);

            for j = 1:181
                o = A(new_array(1,j),1);
                l = A(new_array(1,j),2);
                interaction_matrix1(o,l) = 0;
                interaction_matrix1_index(o,l) = 0;
            end

            protein_gauss = similarity_protein(interaction_matrix1);
            disease_gauss = similarity_disease(interaction_matrix1);
            [protein_integration_similarity,disease_integration_similarity] = ...
                integration_protein_disease_similarity(protein_gauss, disease_gauss);

            predict_score_matrix1 = MDHGIMDA(interaction_matrix1);

            % BNNR 预测的关联得分
            Wdd = protein_integration_similarity;
            Wvv = disease_integration_similarity;
            Wvd = interaction_matrix1';
            [dn,dr] = size(Wvd);

            T = [Wdd, Wvd'; Wvd, Wvv];
            [t1,t2] = size(T);
            trIndex = double(T ~= 0);
            [WW,iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
            M_recovery = WW((t1-dn+1) : t1, 1 : dr);
            predict_score_matrix2 = M_recovery';

            predict_score_matrix3_3 = ensemble_average( ...
                predict_score_matrix1, ...
                predict_score_matrix2, ...
                interaction_matrix1_index, ...
                0.3);

            % 调用第一种方法 RBM
            Sco1 = predict_score_matrix1;
            final_score = Sco1(interaction_matrix == 0);
            for i = 1:181
                q = A(new_array(1,i),1);
                w = A(new_array(1,i),2);
                s_score = Sco1(q,w);
                T = [s_score;final_score];
                index1 = sort(T,'descend');
                index1_1 = find(index1 == s_score);

                if length(index1_1) ~= 1
                    index1_2 = mean(index1_1);
                else
                    index1_2 = index1_1;
                end
                five_g1_1(i) = index1_2;
            end

            % 第二种：INBM 方法
            Sco2 = predict_score_matrix2;
            final_score = Sco2(interaction_matrix == 0);
            for i = 1:181
                q = A(new_array(1,i),1);
                w = A(new_array(1,i),2);
                s_score = Sco2(q,w);
                T = [s_score;final_score];
                index2 = sort(T,'descend');
                index2_1 = find(index2 == s_score);
                if length(index2_1) ~= 1
                    index2_2 = mean(index2_1);
                else
                    index2_2 = index2_1;
                end
                five_g2_1(i) = index2_2;
            end

            five3 = conputing_ensemble_ranking( ...
                predict_score_matrix3_3, ...
                interaction_matrix, ...
                181, ...
                A, ...
                new_array);

            C1_ensemble3 = [C1_ensemble3,five3];
        end
    end

    five_k_ensemble3(k,:) = C1_ensemble3;
    C1_ensemble3 = [];
end


xlswrite(fullfile(result_dir, 'five_fold_ensemble3.xlsx'), five_k_ensemble3);

toc;
fprintf('five_fold结束\n');


interaction_file = fullfile(result_dir, 'Protein_Disease_Associations.xlsx');
adj_file = fullfile(result_dir, 'Protein_Disease_adj.txt');
output3_file = fullfile(result_dir, 'output3.xlsx');

interaction = xlsread(interaction_file);
[n, m] = size(interaction);

adj_pairs = load(adj_file);
[pp, ~] = size(adj_pairs);

files = {
    'five_fold_ensemble3.xlsx'
};

model_names = {
    'ensemble3'
};

mean_auc = zeros(length(files), 1);

for f = 1:length(files)
    file_path = fullfile(result_dir, files{f});
    position_matrix = xlsread(file_path);

    [round_num, ~] = size(position_matrix);
    auc_each_round = nan(round_num, 1);

    for k = 1:round_num
        position = position_matrix(k, :);

        % 如果某些行全是 0，说明这一轮没跑，跳过
        if all(position == 0)
            continue;
        end

        auc_each_round(k) = calc_auc_from_position(position, n, m, pp);
    end

    mean_auc(f) = mean(auc_each_round, 'omitnan');

    fprintf('%s mean AUC = %.6f\n', model_names{f}, mean_auc(f));

    
    auc_output_file = fullfile(result_dir, [model_names{f}, '_auc_each_round.xlsx']);
    xlswrite(auc_output_file, auc_each_round);
end


summary_table = table(model_names, mean_auc, ...
    'VariableNames', {'Model', 'Mean_AUC'});

summary_file = fullfile(result_dir, 'mean_auc_summary.xlsx');
writetable(summary_table, summary_file);

fprintf('AUC 计算完成，结果已保存到：\n%s\n', summary_file);






function overallauc = calc_auc_from_position(position, n, m, pp)

    position = position(:)';

    fold_size = floor(pp / 5);

    max_k = m * n - fold_size * 4;

    denominator = fold_size * 4 * ...
        (m * n - pp + fold_size - 1) + ...
        
        (pp - fold_size * 4) * ...
        (m * n - fold_size * 4 - 1);

    tpr = zeros(max_k, 1);
    fpr = zeros(max_k, 1);

    for kk = 1:max_k
        tp = sum(position <= kk);

        tpr(kk) = tp / pp;

        fp = kk * pp - tp;

        fpr(kk) = fp / denominator;
    end

    overallauc = trapz(fpr, tpr);
end
