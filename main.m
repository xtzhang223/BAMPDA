clear;
clc;

%% 路径设置
result_dir = 'C:\Users\zxt\Desktop\run_MDHGIBNNR\demo';

%% 读取基础数据
A = load(fullfile(result_dir, 'Protein_Disease_adj.txt'));
nd = max(A(:,1));
nv = max(A(:,2));
[pp, qq] = size(A);

interaction_matrix = xlsread(fullfile(result_dir, 'Protein_Disease_Associations.xlsx'));
[row, col] = size(interaction_matrix);

number = 100;

%% BNNR 模型参数
maxiter = 300;
alpha = 1;
beta = 10;
tol1 = 2 * 1e-3;
tol2 = 1 * 1e-5;

array1 = xlsread(fullfile(result_dir, 'array_P-D.xlsx'));

%% 设置权重 w = 0.1 到 0.9
ensemble_weights = 0.1:0.1:0.9;

mean_auc_all = zeros(length(ensemble_weights), 1);
model_names = cell(length(ensemble_weights), 1);

tic;

%% 对每一个权重分别运行
for weight_id = 1:length(ensemble_weights)

    ensemble_weight = ensemble_weights(weight_id);

    weight_tag = sprintf('w%02d', round(ensemble_weight * 10));
    model_name = ['ensemble_', weight_tag];
    model_names{weight_id} = model_name;

    fprintf('\n=====================================\n');
    fprintf('开始运行 ensemble weight = %.1f\n', ensemble_weight);
    fprintf('=====================================\n');

    C1_ensemble3 = [];
    five_k_ensemble3 = zeros(number, pp);

    %% 五折交叉验证
    for k = 95:number

        fprintf('第%d轮\n', k);
        current_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
        disp(current_time);

        array = array1(k,:);

        for circle = 0:4

            if circle < 4

                nu = 181;
                interaction_matrix1 = interaction_matrix;
                interaction_matrix1_index = logical(zeros(row, col) ~= 0);

                for i = 1:pp
                    interaction_matrix1_index(A(i,1), A(i,2)) = 1;
                end

                new_array = array(1, 1 + circle * 181 : (circle + 1) * 181);

                for j = 1:181
                    o = A(new_array(1,j), 1);
                    l = A(new_array(1,j), 2);
                    interaction_matrix1(o,l) = 0;
                    interaction_matrix1_index(o,l) = 0;
                end

                protein_gauss = xlsread(fullfile(result_dir, 'Protein_Gaussian_Similarity_Matrix.xlsx'));
                disease_gauss = xlsread(fullfile(result_dir, 'Disease_Gaussian_Similarity_Matrix.xlsx'));

                [protein_integration_similarity, disease_integration_similarity] = ...
                    integration_protein_disease_similarity(protein_gauss, disease_gauss);

                %% MDHGIMDA 预测得分
                predict_score_matrix1 = MDHGIMDA(interaction_matrix1);

                %% BNNR 预测得分
                Wdd = protein_integration_similarity;
                Wvv = disease_integration_similarity;
                Wvd = interaction_matrix1';

                [dn, dr] = size(Wvd);

                T = [Wdd, Wvd'; Wvd, Wvv];
                [t1, t2] = size(T);
                trIndex = double(T ~= 0);

                [WW, iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);

                M_recovery = WW((t1 - dn + 1) : t1, 1 : dr);
                predict_score_matrix2 = M_recovery';

                %% ensemble 预测，权重为 ensemble_weight
                predict_score_matrix3_3 = ensemble_average( ...
                    predict_score_matrix1, ...
                    predict_score_matrix2, ...
                    interaction_matrix1_index, ...
                    ensemble_weight);

                %% 计算 ensemble 排名
                five3 = conputing_ensemble_ranking( ...
                    predict_score_matrix3_3, ...
                    interaction_matrix, ...
                    181, ...
                    A, ...
                    new_array);

                C1_ensemble3 = [C1_ensemble3, five3];

            end

            if circle == 4

                nu = 181;
                interaction_matrix1 = interaction_matrix;
                interaction_matrix1_index = logical(zeros(row, col) ~= 0);

                for i = 1:pp
                    interaction_matrix1_index(A(i,1), A(i,2)) = 1;
                end

                new_array = array(1, 1 + circle * 181 : pp);

                for j = 1:181
                    o = A(new_array(1,j), 1);
                    l = A(new_array(1,j), 2);
                    interaction_matrix1(o,l) = 0;
                    interaction_matrix1_index(o,l) = 0;
                end

                protein_gauss = similarity_protein(interaction_matrix1);
                disease_gauss = similarity_disease(interaction_matrix1);

                [protein_integration_similarity, disease_integration_similarity] = ...
                    integration_protein_disease_similarity(protein_gauss, disease_gauss);

                %% MDHGIMDA 预测得分
                predict_score_matrix1 = MDHGIMDA(interaction_matrix1);

                %% BNNR 预测得分
                Wdd = protein_integration_similarity;
                Wvv = disease_integration_similarity;
                Wvd = interaction_matrix1';

                [dn, dr] = size(Wvd);

                T = [Wdd, Wvd'; Wvd, Wvv];
                [t1, t2] = size(T);
                trIndex = double(T ~= 0);

                [WW, iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);

                M_recovery = WW((t1 - dn + 1) : t1, 1 : dr);
                predict_score_matrix2 = M_recovery';

                %% ensemble 预测，权重为 ensemble_weight
                predict_score_matrix3_3 = ensemble_average( ...
                    predict_score_matrix1, ...
                    predict_score_matrix2, ...
                    interaction_matrix1_index, ...
                    ensemble_weight);

                %% 计算 ensemble 排名
                five3 = conputing_ensemble_ranking( ...
                    predict_score_matrix3_3, ...
                    interaction_matrix, ...
                    181, ...
                    A, ...
                    new_array);

                C1_ensemble3 = [C1_ensemble3, five3];

            end

        end

        five_k_ensemble3(k,:) = C1_ensemble3;
        C1_ensemble3 = [];

    end

    %% 保存当前权重的五折排名结果
    five_fold_file = fullfile(result_dir, ...
        ['five_fold_ensemble_', weight_tag, '.xlsx']);

    xlswrite(five_fold_file, five_k_ensemble3);

    fprintf('当前权重 %.1f 的五折结果已保存到：\n%s\n', ...
        ensemble_weight, five_fold_file);

    %% 计算当前权重的 AUC
    interaction = xlsread(fullfile(result_dir, 'Protein_Disease_Associations.xlsx'));
    [n, m] = size(interaction);

    position_matrix = xlsread(five_fold_file);
    [round_num, ~] = size(position_matrix);

    auc_each_round = nan(round_num, 1);

    for k = 1:round_num

        position = position_matrix(k, :);

        if all(position == 0)
            continue;
        end

        auc_each_round(k) = calc_auc_from_position(position, n, m, pp);

    end

    mean_auc = mean(auc_each_round, 'omitnan');
    mean_auc_all(weight_id) = mean_auc;

    fprintf('w = %.1f, mean AUC = %.6f\n', ensemble_weight, mean_auc);

    %% 保存当前权重每一轮的 AUC
    auc_output_file = fullfile(result_dir, ...
        [model_name, '_auc_each_round.xlsx']);

    xlswrite(auc_output_file, auc_each_round);

    fprintf('当前权重 %.1f 的每轮 AUC 已保存到：\n%s\n', ...
        ensemble_weight, auc_output_file);

end

%% 汇总保存所有权重的 mean AUC
summary_table = table(ensemble_weights', model_names, mean_auc_all, ...
    'VariableNames', {'Weight', 'Model', 'Mean_AUC'});

summary_file = fullfile(result_dir, 'mean_auc_summary_all_weights.xlsx');
writetable(summary_table, summary_file);

fprintf('\n所有权重运行完成，AUC 汇总结果已保存到：\n%s\n', summary_file);

toc;
fprintf('全部实验结束。\n');


%% ====== 本地函数：根据一轮排名计算 AUC ======
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