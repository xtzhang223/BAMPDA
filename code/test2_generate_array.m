num_groups = 100;  % 需要生成 100 组
num_range = 905;  % 数值范围 1-1000

all_sequences = zeros(num_groups, num_range); % 预分配存储空间，提高效率

for i = 1:num_groups
    all_sequences(i, :) = randperm(num_range); % 生成 1 到 1000 的随机排列
end

% 显示部分数据
disp(all_sequences(1:5, 1:10)); % 仅显示前 5 组的前 10 个数，避免输出过长
