function [t, x, v] = solve_liquid_bridge(beta, Vol_uL, T_grad, Phi_deg, mu_coeff)
% SOLVE_LIQUID_BRIDGE 求解液桥热毛细迁移的动力学方程
%
% 输入参数:
%   beta      : 接触角滞后系数 (Hysteresis coefficient)
%   Vol_uL    : 液桥体积 (微升 uL)
%   T_grad    : 温度梯度 (K/m)
%   Phi_deg   : 倾斜角度 (度)
%   mu_coeff  : 粘温系数 (1/K)
%
% 输出参数:
%   t : 时间向量 (s)
%   x : 位移向量 (m)
%   v : 速度向量 (m/s)

    %% 1. 固定物理参数定义 (Fixed Parameters)
    mu0     = 0.105;          % 热端初始粘度 [Pa.s] (100 mPa.s)
    rho     = 963;          % 密度 [kg/m^3]
    gamma_0 = 0.021;        % 热端表面张力 [N/m]
    gamma_T = 0.046e-3;     % 表面张力温度系数 [N/m.K]
    
    H_gap   = 1.25e-3;      % 平板间距 2h [m]
    h       = H_gap / 2;    % 半高 h [m]
    g       = 9.81;         % 重力加速度 [m/s^2]
    
    % 接触角 (固定为 11.5 度)
    theta_a = 11.5 * (pi/180); 

    %% 2. 几何与单位换算
    Vol = Vol_uL * 1e-9;         % 换算为 [m^3]
    Phi = Phi_deg * (pi/180);    % 换算为 [rad]
    
    % 计算特征直径 D
    D = sqrt(2 * Vol / (pi * h));
    
    %% 3. 计算常数系数 (Pre-calculate Constants)
    
    % Term 1: Marangoni Drive (热毛细驱动项)
    C_drive = (4 * gamma_T * T_grad * cos(theta_a)) / (pi * rho * h);
    
    % Term 2: Hysteresis Factor (滞后阻力项系数)
    C_hyst_factor = (4 * beta * cos(theta_a)) / (pi * rho * h * D);
    
    % Term 3: Drag Base Factor (粘性阻力项系数 - 理论值12，可根据需要修正)
    C_drag_base = 4 / (pi * rho * h^2);
    
    % Term 4: Gravity (重力项)
    Acc_grav = g * sin(Phi);
    
    %% 4. 配置 ODE 求解器
    y0 = [0; 0];       % 初始状态 [位置; 速度]
    tspan = [0 80];    % 仿真时间 [s]
    
    % 设置精度
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    
    % 定义匿名函数，传递参数
    ode_system = @(t, y) dynamics_func(t, y, mu0, mu_coeff, T_grad, gamma_0, gamma_T, ...
                                       C_drive, C_hyst_factor, C_drag_base, Acc_grav);
    
    % 调用 ode45
    [t, y] = ode45(ode_system, tspan, y0, options);
    
    % 提取结果
    x = y(:, 1);
    v = y(:, 2);

end

%% 5. 局部动力学方程 (Local Dynamics Function)
function dydt = dynamics_func(~, y, mu0, mu_coeff, T_grad, gamma_0, gamma_T, ...
                              C_drive, C_hyst_factor, C_drag_base, Acc_grav)
    % 解包状态
    x_curr = y(1);
    v_curr = y(2);
    
    % A. 更新随位置变化的物性
    % 粘度随位置指数变化
    mu_x = mu0 * exp(mu_coeff * T_grad * x_curr);
    % 表面张力随位置线性变化
    gamma_x = gamma_0 + gamma_T * T_grad * x_curr;
    
    % B. 计算加速度
    % 静态加速度 = 驱动力 - 滞后阻力 - 重力
    acc_static = C_drive - (C_hyst_factor * gamma_x) - Acc_grav;
    
    % 粘性阻力加速度 = - C_drag * mu(x) * v
    acc_drag = (C_drag_base * mu_x) * v_curr;
    
    % C. 总加速度 (dv/dt)
    dvdt = acc_static - acc_drag;
    
    % 物理锁定逻辑 (防止反向滑动): 
    % 如果速度<=0 且 驱动力不足以克服静态阻力，则保持静止
    if v_curr <= 0 && acc_static <= 0
        dvdt = 0;
        v_curr = 0;
    end
    
    % D. 打包输出 [dx/dt; dv/dt]
    dydt = [v_curr; dvdt];
end