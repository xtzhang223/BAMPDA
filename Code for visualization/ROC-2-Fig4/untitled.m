% 设置文件夹路径
folder_path = 'C:\Users\zxt\Desktop\run_MDHGIBNNR\ROC-2';
desktop_path = 'C:\Users\zxt\Desktop\run_MDHGIBNNR';

% 获取所有 CSV 文件
files = dir(fullfile(folder_path, '*.csv'));

% 颜色列表
colors = {'b', 'g', 'r', 'c', 'm', 'y'};

figure; hold on;

% 设置全局字体（Times New Roman + 字号）
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
set(gcf, 'Color', 'w');

% 依次读取 CSV 文件并绘制 ROC 曲线
for i = 1:length(files)
    file_path = fullfile(folder_path, files(i).name);
    data = readmatrix(file_path);
    fpr = data(:, 1);
    tpr = data(:, 2);

    plot(fpr, tpr, 'Color', colors{mod(i-1, length(colors)) + 1}, ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('ROC %d (%s)', i, files(i).name));
end

% 添加图例和标签（Times New Roman）
xlabel('False Positive Rate (FPR)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('True Positive Rate (TPR)', 'FontSize', 20, 'FontName', 'Times New Roman');
title('ROC curves based on 5-CV experiments', 'FontSize', 22, 'FontName', 'Times New Roman');
legend('Location', 'southeast', 'FontSize', 18, 'FontName', 'Times New Roman');

grid on;
hold off;

% 设置图形大小
set(gcf, 'Position', [100, 100, 800, 600]);

% 保存为 JPG 格式
output_filename = fullfile(desktop_path, 'ROC_curves_400dpi.jpg');
print(output_filename, '-djpeg', '-r400');

disp(['图像已保存到桌面: ' output_filename]);
