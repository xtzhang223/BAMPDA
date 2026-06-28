function [overallauc, precision, recall, f1_score] = calculate_metrics_ensemble3_original()

clear;
clc;

%% 路径设置
data_dir = 'C:\Users\zxt\Desktop\run_MDHGIBNNR\demo';

position_file = fullfile(data_dir, 'five_fold_ensemble3.xlsx');
true_matrix_file = fullfile(data_dir, 'Protein_Disease_Associations.xlsx');
predicted_matrix_file = fullfile(data_dir, 'output1.xlsx');

adj_xlsx_file = fullfile(data_dir, 'Protein_Disease_adj .xlsx');
adj_txt_file = fullfile(data_dir, 'Protein_Disease_adj.txt');

%% =========================
%  1. AUC 计算：保持原始 positiontooverallaucfold 逻辑
%  =========================

position1 = xlsread(position_file);

% 原始代码中 for i=1:100 最后实际只保留第100行 position
for i = 1:100
    position = position1(i, :);

    interaction = xlsread(true_matrix_file);
    [n, m] = size(interaction);

    if exist(adj_txt_file, 'file')
        sID = textread(adj_txt_file);
    else
        sID = xlsread(adj_xlsx_file);
    end

    [pp, qq] = size(sID);
end

for k = 1:m * n - floor(pp / 5) * 4
    tp = 0;

    for t = 1:pp
        if position(1, t) <= k
            tp = tp + 1;
        end
    end

    tpr(1, k) = tp / pp;
    fp = k * pp - tp;

    fpr(1, k) = fp / ...
        (floor(pp / 5) * 4 * (m * n - pp + floor(pp / 5) - 1) + ...
        (pp - floor(pp / 5) * 4) * (m * n - floor(pp / 5) * 4 - 1));
end

overallauc = trapz(fpr, tpr);

figure;
plot(fpr, tpr, 'LineWidth', 2);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(['ROC curve, AUC = ', num2str(overallauc, '%.4f')]);
grid on;

%% =========================
%  2. Precision / Recall / F1-score
%  保持原始 0/1 矩阵计算公式
%  =========================

A = readmatrix(true_matrix_file);
B = readmatrix(predicted_matrix_file);

if ~isequal(size(A), size(B))
    error('矩阵维度不一致，请检查数据！');
end

A = logical(A);
B = logical(B);

TP = sum(A & B, 'all');
FP = sum(~A & B, 'all');
FN = sum(A & ~B, 'all');

if (TP + FP) == 0
    precision = 0;
else
    precision = TP / (TP + FP);
end

if (TP + FN) == 0
    recall = NaN;
else
    recall = TP / (TP + FN);
end

if precision + recall == 0
    f1_score = 0;
else
    f1_score = 2 * precision * recall / (precision + recall);
end

%% 输出结果
fprintf('AUC:       %.4f\n', overallauc);
fprintf('Recall:    %.4f\n', recall);
fprintf('Precision: %.4f\n', precision);
fprintf('F1-score:  %.4f\n', f1_score);

%% 保存结果
metric_result = {
    'Metric', 'Value';
    'AUC', overallauc;
    'Recall', recall;
    'Precision', precision;
    'F1-score', f1_score
};

xlswrite(fullfile(data_dir, 'metrics_ensemble3_original.xlsx'), metric_result);

end