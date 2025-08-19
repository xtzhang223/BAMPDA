
vv = xlsread('virus sequence similarity.xlsx'); % 示例稠密矩阵
c = vv;
% 将矩阵中的元素按绝对值排序
sorted_values_vv = sort(abs(vv(:)), 'descend');



% 假设 ss 是一个稠密的相似性矩阵
ss = xlsread('Protein_Sequence_Similarity_Matrix.xlsx'); % 示例稠密矩阵
a = ss;
threshold_percentile = 70;  % 设置90%的百分位数作为阈值

% 将矩阵中的元素按绝对值排序
sorted_values_ss = sort(abs(ss(:)), 'descend');

% 选择阈值
threshold_value = sorted_values_ss(round(threshold_percentile / 100 * length(sorted_values_ss)));

% 将小于阈值的元素设为零
ss(ss < threshold_value) = 0;

dd = xlsread('Disease_Semantic_Similarity_Matrix.xlsx'); % 示例稠密矩阵
b = dd;


% 将矩阵中的元素按绝对值排序
sorted_values_dd = sort(abs(dd(:)), 'descend');

% 选择阈值
threshold_value = sorted_values_dd(round(threshold_percentile / 100 * length(sorted_values_dd)));

% 将小于阈值的元素设为零
dd(dd < threshold_value) = 0;

xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\代码\pp_decline.xlsx',ss);
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\代码\disease_decline.xlsx',dd);



% 输出稀疏矩阵
disp(ss);
