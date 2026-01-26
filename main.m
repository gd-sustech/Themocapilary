% =========================================================================
% MATLAB Code: Dynamics of Liquid Bridge (Theta_a Model) using ode45
% Model: Hot (Bottom) -> Cold (Top) [Uphill]
% Equation Reference: "Dynamics of Thermocapillary Migration with Positive Hysteresis"
% Method: Runge-Kutta (4,5) Formula (ode45)
% =========================================================================

clear; clc; close all;

%% 1. Parameter Definition

% --- Liquid Properties (Silicone Oil) ---
mu0     = 0.1;          % Dynamic Viscosity at Hot End [Pa.s]
rho     = 963;          % Density [kg/m^3]
gamma_0 = 0.021;         % Surface Tension at Hot End [N/m]
gamma_T = 0.046e-3;     % Surface Tension Temp Coeff [N/m.K]
mu_coeff= 0;         % Viscosity-Temp Coefficient [1/K]

% --- Contact Angle & Hysteresis ---
theta_a = 11.5 * (pi/180);% Advancing Contact Angle [rad] (Cold end)
beta    = 0.0;            % Hysteresis coefficient (Resistance)

% --- Geometry & Conditions ---
Vol_uL  = 20;           % Volume [uL]
Vol     = Vol_uL * 1e-9;% [m^3]
H_gap   = 1.25e-3;      % Gap height 2h [m]
h       = H_gap / 2;    % Half-height [m]
T_grad  = 3 * 1e3;      % Thermal Gradient [K/m]

% --- Inclination (Gravity) ---
Phi_deg = 0;            % Tilt Angle [deg] (Positive = Uphill)
Phi     = Phi_deg * (pi/180);
g       = 9.81;         % Gravity [m/s^2]

% Calculate Characteristic Diameter D
D = sqrt(2 * Vol / (pi * h));

fprintf('--- Simulation Configuration (ode45) ---\n');
fprintf('Model: Reference Angle theta_a = %.1f deg\n', theta_a*(180/pi));
fprintf('Hysteresis beta: %.3f\n', beta);
fprintf('Tilt Angle: %.1f deg\n', Phi_deg);
fprintf('Diameter D: %.2f mm\n', D*1e3);

%% 2. Solver Setup

% --- Pre-calculate Constants (Pass these to the function) ---
% Term 1: Marangoni Drive
C_drive = (4 * gamma_T * T_grad * cos(theta_a)) / (pi * rho * h);

% Term 2: Hysteresis Factor
C_hyst_factor = (4 * beta * cos(theta_a)) / (pi * rho * h * D);

% Term 3: Drag Base Factor
C_drag_base = 12 / (pi * rho * h^2);

% Term 4: Gravity
Acc_grav = g * sin(Phi);

% Initial State Vector: [Position; Velocity]
y0 = [0; 0]; 

% Time Span
tspan = [0 80];

% Pack parameters into a struct or pass directly to nested function
% Using nested function approach for simplicity in single file
ode_system = @(t, y) dynamics_func(t, y, mu0, mu_coeff, T_grad, gamma_0, gamma_T, ...
                                   C_drive, C_hyst_factor, C_drag_base, Acc_grav);

%% 3. Solve using ode45
% Options: RelTol and AbsTol control accuracy
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

[t, y] = ode45(ode_system, tspan, y0, options);

% Extract results
x = y(:, 1);
v = y(:, 2);

%% 4. Visualization
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

% Plot 1: Displacement (x-t)
subplot(2,1,1);
plot(t, x * 1e3, 'b-', 'LineWidth', 2);
xlabel('Time t (s)', 'FontSize', 12);
ylabel('Position x (mm)', 'FontSize', 12);
title(['Migration Distance (\theta_a=' num2str(theta_a*180/pi) '^{\circ}, \beta=' num2str(beta) ')'], 'FontSize', 14);
grid on;

% Plot 2: Velocity (v-t)
subplot(2,1,2);
plot(t, v * 1e3, 'r-', 'LineWidth', 2);
xlabel('Time t (s)', 'FontSize', 12);
ylabel('Velocity v (mm/s)', 'FontSize', 12);
title('Migration Velocity', 'FontSize', 14);
grid on;

% Results Output
fprintf('\n--- Simulation Results ---\n');
fprintf('Final Position: %.2f mm\n', x(end)*1e3);
fprintf('Max Velocity:   %.2f mm/s\n', max(v)*1e3);

%% 5. Local Function Definition
function dydt = dynamics_func(~, y, mu0, mu_coeff, T_grad, gamma_0, gamma_T, ...
                              C_drive, C_hyst_factor, C_drag_base, Acc_grav)
    % Unpack state
    x_curr = y(1);
    v_curr = y(2);
    
    % A. Update Position-Dependent Properties
    % Viscosity increases with x
    mu_x = mu0 * exp(mu_coeff * T_grad * x_curr);
    % Surface Tension increases with x
    gamma_x = gamma_0 + gamma_T * T_grad * x_curr;
    
    % B. Calculate Accelerations
    % Static Acceleration (Drive - Hysteresis - Gravity)
    acc_static = C_drive - (C_hyst_factor * gamma_x) - Acc_grav;
    
    % Viscous Drag Acceleration ( - C_drag * mu * v )
    acc_drag = (C_drag_base * mu_x) * v_curr;
    
    % C. Total Acceleration (dv/dt)
    dvdt = acc_static - acc_drag;
    
    % D. Pack output (dx/dt = v; dv/dt = a)
    dydt = [v_curr; dvdt];
end