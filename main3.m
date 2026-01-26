clear; clc; close all;


excel_file = 'experiment_data.xlsx'; % 请确保文件名一致

try
    % 读取整个矩阵 (假设有8列数据: t1, x1, t2, x2, t3, x3, t4, x4)
    raw_exp_data = readmatrix(excel_file);
    fprintf('成功读取 Excel 文件: %s\n', excel_file);
catch
    error('未找到文件 "%s"。请创建该文件并填入实验数据，或修改代码中的文件名。', excel_file);
end


% {beta, Vol_uL, T_grad, Phi_deg, mu_coeff, 'Legend Name'}
cases = {
    {0.013,  20, 5900, 0, 0.00, 'Case 1: 20uL, 5.9K/mm'},
    {0.013,  20, 4300, 0, 0.0, 'Case 2: 20uL, 4.3K/mm'},
    {0.013,  20, 3000, 0, 0.0, 'Case 3: 20uL, 3.0K/mm'},
    {0.01,  40, 3000, 0, 0.0, 'Case 4: 40uL, 3.0K/mm'}
};

styles = {'-', '-', '-', '-'}; 
colors = {'m', 'r', 'k', 'b'}; 
data_log = struct(); 


figure('Color', 'w', 'Position', [100, 100, 900, 600]);
hold on;

h_sim = []; h_exp = []; % 用于存储图例句柄

for i = 1:length(cases)
    % --- A. 理论模拟 ---
    [beta, V, G, Phi, mu_c, leg] = cases{i}{:};
    
    [t, x, ~] = solve_liquid_bridge(beta, V, G, Phi, mu_c);
    x_mm = x * 1e3;
    
    % 保存模拟数据
    data_log(i).case_name = leg;
    data_log(i).sim_t = t;
    data_log(i).sim_x_mm = x_mm;
    
    % 绘制模拟曲线
    p_sim = plot(t, x_mm, 'LineStyle', styles{i}, 'Color', colors{i}, ...
        'LineWidth', 2, 'DisplayName', [leg ' (Sim)']);
    h_sim = [h_sim, p_sim];
    
    % --- B. 实验数据处理与绘制 ---
    % 计算对应的列索引: Case 1 -> 1,2列; Case 2 -> 3,4列 ...
    col_idx_t = 2*i - 1;
    col_idx_x = 2*i;
    
    if size(raw_exp_data, 2) >= col_idx_x
        % 提取数据列
        t_e_raw = raw_exp_data(:, col_idx_t);
        x_e_raw = raw_exp_data(:, col_idx_x);
        
        % 去除 NaN (防止不同组数据长度不一致带来的空值干扰)
        valid_idx = ~isnan(t_e_raw) & ~isnan(x_e_raw);
        t_e = t_e_raw(valid_idx);
        x_e = x_e_raw(valid_idx);
        
        % 保存实验数据到结构体
        data_log(i).exp_t = t_e;
        data_log(i).exp_x_mm = x_e;
        
        % 绘制实验点 (*)
        p_exp = plot(t_e, x_e, '*', 'Color', colors{i}, ...
            'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', [leg ' (Exp)']);
        h_exp = [h_exp, p_exp];
    else
        warning('Excel 中缺少 Case %d 的数据列 (需要第 %d-%d 列)', i, col_idx_t, col_idx_x);
    end
end


xlabel('Time t (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Position x (mm)', 'FontSize', 12, 'FontWeight', 'bold');
title('Comparison: Simulation (-) vs Experiment (*)', 'FontSize', 14);

grid on; box on;
xlim([0, 80]); 
% 自动调整 Y 轴上限
all_y = [cellfun(@(x) max(x), {data_log.sim_x_mm})];
if isfield(data_log, 'exp_x_mm')
    for k=1:length(data_log)
        if ~isempty(data_log(k).exp_x_mm)
            all_y = [all_y, max(data_log(k).exp_x_mm)];
        end
    end
end
ylim([0, max(all_y) * 1.15]);

% 图例
legend([h_sim, h_exp], 'Location', 'best', 'FontSize', 9, 'NumColumns', 2);

hold off;

% 保存结果
save('comparison_results.mat', 'data_log');
fprintf('所有模拟与实验数据已保存至 "comparison_results.mat"。\n');