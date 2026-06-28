clear;
clc;

A = xlsread('C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\Protein_Disease_adj.xlsx');  
A = A(:,1:2);  % 仅保留索引信息，确保维度一致

nd = max(A(:,1));       % protein 数
nv = max(A(:,2));       % disease 数
[pp, qq] = size(A);     % 已知正样本数

% 读取原始互动矩阵：1表示已有关联，0表示无关联
interaction_matrix = xlsread('C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\Protein_Disease_Associations.xlsx');
[row, col] = size(interaction_matrix);

number = 100;           % 外层循环次数

% BNNR 模型参数
maxiter = 300;
alpha   = 1;
beta    = 10;
tol1    = 2*1e-3;
tol2    = 2*1e-5;

% 读取用于交叉验证的随机排列数组
array1 = xlsread('C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\array_P-D.xlsx');  
% array1 为 number x pp 矩阵，每行为正样本 A 的一个随机排列

[xg, yg] = ndgrid(1:nd, 1:nv);
all_pairs = [xg(:), yg(:)];
global_neg_pool = setdiff(all_pairs, A, 'rows');  % 全局负样本池

tic;

for k = 95:number
    fprintf('第%d轮\n', k);
    current_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
    disp(current_time);
    
    % 从随机排列数组中读取当前轮的排列
    array = array1(k, :);
    
    for circle = 0:4
        % 每折隐藏的正样本数
        nu = 181;

        % 复制原始互动矩阵作为当前折的训练矩阵
        interaction_matrix1 = interaction_matrix;
        interaction_matrix1_index = false(row, col);

        % 构造正样本索引矩阵：初始将所有正样本置 true
        for i = 1:pp
            interaction_matrix1_index(A(i,1), A(i,2)) = true;
        end

        % 选取每折测试集的正样本索引
        new_array = array(1, 1+circle*nu : (circle+1)*nu);
        test_positive = A(new_array, :);

        % 隐藏测试集的正样本
        for j = 1:nu
            o = A(new_array(1,j), 1);
            l = A(new_array(1,j), 2);

            interaction_matrix1(o, l) = 0;
            interaction_matrix1_index(o, l) = false;
        end

        % 构造训练集的正样本
        % 从剩余正样本中随机选取 724 个样本
        remaining_positive_indices = setdiff(1:pp, new_array);

        perm_pos = randperm(length(remaining_positive_indices), 724);
        train_positive_indices = remaining_positive_indices(perm_pos);

        train_positive = A(train_positive_indices, :);

        % 从全局负样本池中随机采样 724 个负样本
        neg_perm = randperm(size(global_neg_pool, 1), 724);
        train_negative = global_neg_pool(neg_perm, :);

        % 训练集的正负样本合并
        train_set = [train_positive, ones(size(train_positive, 1), 1);
                     train_negative, zeros(size(train_negative, 1), 1)];

        % 测试集仍按原逻辑构造
        test_negative_indices = randperm(size(global_neg_pool, 1), nu);
        test_negative = global_neg_pool(test_negative_indices, :);

        balanced_test_set = [test_positive, ones(size(test_positive, 1), 1);
                             test_negative, zeros(size(test_negative, 1), 1)];

        % 将采样的负样本也显式标记为 0
        for i = 1:size(train_negative, 1)
            r = train_negative(i, 1);
            c = train_negative(i, 2);

            interaction_matrix1(r, c) = 0;
            interaction_matrix1_index(r, c) = false;
        end
    end

    % ---------------------------------------------------------

    % 加载相似性矩阵
    protein_gauss = xlsread('C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\Protein_Gaussian_Similarity_Matrix.xlsx');
    disease_gauss = xlsread('C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\Disease_Gaussian_Similarity_Matrix.xlsx');

    [protein_integration_similarity, disease_integration_similarity] = ...
        integration_protein_disease_similarity(protein_gauss, disease_gauss);
    
    % MDHGIMDA 模型预测
    predict_score_matrix1 = MDHGIMDA(interaction_matrix1);
    
    % BNNR 模型预测
    Wdd = protein_integration_similarity;
    Wvv = disease_integration_similarity;
    Wvd = interaction_matrix1';

    [dn, dr] = size(Wvd);

    T = [Wdd, Wvd'; Wvd, Wvv];
    [t1, ~] = size(T);
    trIndex = double(T ~= 0);

    [WW, iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);

    M_recovery = WW((t1-dn+1):t1, 1:dr);
    predict_score_matrix2 = M_recovery';
    
    % 集成 MDHGIMDA 和 BNNR
    predict_score_matrix3_3 = ensemble_average( ...
        predict_score_matrix1, ...
        predict_score_matrix2, ...
        interaction_matrix1_index, ...
        0.3);
end

%% 保存连续预测得分矩阵
xlswrite( ...
    'C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\Predict_Score_Matrix3_3.xlsx', ...
    predict_score_matrix3_3);

%% 固定阈值生成0/1预测矩阵
threshold = 0.393178;

% 大于阈值为1，否则为0
binary_matrix = double(predict_score_matrix3_3 > threshold);

% 统计最终1的数量
num_positive = sum(binary_matrix(:) == 1);

%% 保存0/1预测矩阵
xlswrite( ...
    'C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\output3.xlsx', ...
    binary_matrix);

%% 保存阈值
xlswrite( ...
    'C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\threshold_fixed.xlsx', ...
    threshold);

%% 保存阈值和1的数量表格
threshold_table = table(threshold, num_positive, ...
    'VariableNames', {'Threshold', 'Number_of_Ones'});

writetable( ...
    threshold_table, ...
    'C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\threshold_fixed_table.xlsx');

%% 计算 Precision、Recall 和 F1-score
true_matrix = interaction_matrix;

% 检查维度是否一致
if ~isequal(size(true_matrix), size(binary_matrix))
    error('真实矩阵和预测矩阵维度不一致，请检查数据！');
end

true_binary = logical(true_matrix);
pred_binary = logical(binary_matrix);

TP = sum(true_binary & pred_binary, 'all');
FP = sum(~true_binary & pred_binary, 'all');
FN = sum(true_binary & ~pred_binary, 'all');
TN = sum(~true_binary & ~pred_binary, 'all');

if (TP + FP) == 0
    precision = 0;
else
    precision = TP / (TP + FP);
end

if (TP + FN) == 0
    recall = 0;
else
    recall = TP / (TP + FN);
end

if (precision + recall) == 0
    f1_score = 0;
else
    f1_score = 2 * precision * recall / (precision + recall);
end

%% 输出指标
fprintf('Predict_Score_Matrix3_3.xlsx 已保存。\n');
fprintf('output3.xlsx 已保存。\n');
fprintf('Precision = %.6f\n', precision);
fprintf('Recall    = %.6f\n', recall);
fprintf('F1-score  = %.6f\n', f1_score);

%% 保存 Precision、Recall、F1-score
metric_table = table(threshold, TP, FP, FN, TN, precision, recall, f1_score, ...
    'VariableNames', {'Threshold', 'TP', 'FP', 'FN', 'TN', 'Precision', 'Recall', 'F1_score'});

writetable( ...
    metric_table, ...
    'C:\Users\zxt\Desktop\run_MDHGIBNNR\demo\precision_recall_f1_summary.xlsx');

toc;
fprintf('five_fold结束\n');