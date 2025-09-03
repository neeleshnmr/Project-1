clear; clc; close all;

comm_params = struct();

comm_params.A_gs = 0.5562;      
comm_params.A_sc = 0.4903;      
comm_params.A_ct = 0.4903 * 0.1;
comm_params.A_wf = 0.4903 * 0.1; 


comm_params.delta_gs = 0.003;   
comm_params.delta_sc = 0.0003;  
comm_params.delta_ca = 0.0015;  
comm_params.delta_ct = 0.01;    
comm_params.delta_wf = 0.01;    


comm_params.rho_gs = 2500;      
comm_params.rho_sc = 2330;      
comm_params.rho_ct = 8960;      
comm_params.rho_wf = 1000;      


comm_params.cp_gs = 750;        
comm_params.cp_sc = 700;        
comm_params.cp_ct = 385;        
comm_params.cp_wf = 4180;       


comm_params.m_dot_wf = 0.02;    % Water mass flow rate [kg/s]
comm_params.T_w_in = 30;        % Inlet water temperature [°C]

% Time-varying environmental data (digitized approximations from Fig. 9)
comm_params.time_data = [0, 3600, 7200, 10800, 14400, 18000, 21600, 25200];  % Seconds
comm_params.G_data = [200, 600, 900, 1100, 1000, 800, 500, 300];  % W/m²
comm_params.T_am_data = [25, 28, 30, 32, 31, 30, 29, 28];  % °C


comm_params.alpha_gs = 0.05;    % Glass absorptance
comm_params.alpha_sc = 0.9;     % Solar cells absorptance
comm_params.E_el = 30;          % Electrical power [W] (simplified)

% Parameters for conv coeff calcs
comm_params.epsilon_g = 0.9;    % Glass emissivity
comm_params.sigma = 5.67e-8;    % Stefan-Boltzmann constant [W/m²·K⁴]
comm_params.Vw = 1;             % Wind speed [m/s] for h_cv_gs_am
comm_params.K_g = 0.96;         % Thermal conductivity glass [W/m·K]
comm_params.K_sc = 148;         % Thermal conductivity solar cells (silicon) [W/m·K]
comm_params.K_ct = 385;         % Thermal conductivity copper tube [W/m·K]


% Copper Absorber
params_cu = comm_params;
params_cu.A_ca = 0.4903;          % Copper absorber area [m²]
params_cu.rho_ca = 8960;          % Copper density [kg/m³]
params_cu.cp_ca = 385;            % Copper specific heat [J/kg·°C]
params_cu.K_ca = 385;             % Copper thermal conductivity [W/m·K]
params_cu.material = 'Copper';

% Aluminum Absorber
params_al = comm_params;
params_al.A_ca = 0.4903;          % Aluminum absorber area [m²] (assume same)
params_al.rho_ca = 2700;          % Aluminum density [kg/m³]
params_al.cp_ca = 900;            % Aluminum specific heat [J/kg·°C]
params_al.K_ca = 237;             % Aluminum thermal conductivity [W/m·K]
params_al.material = 'Aluminum';


function dTdt = heatBalanceODE(t, T, params)
    % Unpack state vector (5 variables)
    Tgs = T(1); Tsc = T(2); Tca = T(3); Tct = T(4); Twf = T(5);
    
    % Interpolate time-varying inputs
    G = interp1(params.time_data, params.G_data, t, 'linear', 'extrap');
    T_am = interp1(params.time_data, params.T_am_data, t, 'linear', 'extrap');
    T_s = 0.0552 * T_am^1.5;  % Sky temp
    
    % dynamic heat transfer coefficients
    h_rd_gs_s = params.epsilon_g * params.sigma * (Tgs^2 + T_s^2) * (Tgs + T_s);
    h_cv_gs_am = 2.8 + 3 * params.Vw;
    h_cd_gs_sc = 1 / (params.delta_gs / params.K_g + params.delta_sc / params.K_sc);
    h_cd_sc_ca = 1 / (params.delta_sc / params.K_sc + params.delta_ca / params.K_ca);
    h_cd_ca_ct = 1 / (params.delta_ca / params.K_ca + params.delta_ct / params.K_ct);
    h_cv_ct_wf = 100;  % Simplified; can add dynamic Nu for water
    
    % absorbed energies
    Egs = params.alpha_gs * G * params.A_gs;  % Eq. 6
    Esc = params.alpha_sc * G * params.A_sc;  % Eq. 7
    
    % Assume outlet temp approx. equal to fluid temp
    T_w_out = Twf;
    
    % Adapted Equations
    dTgs_dt = (Egs - h_rd_gs_s * params.A_gs * (Tgs - T_s) ...
               - h_cv_gs_am * params.A_gs * (Tgs - T_am) ...
               - h_cd_gs_sc * params.A_gs * (Tgs - Tsc)) ...
               / (params.rho_gs * params.A_gs * params.delta_gs * params.cp_gs);
    
    dTsc_dt = (Esc + h_cd_gs_sc * params.A_sc * (Tgs - Tsc) ...
               - h_cd_sc_ca * params.A_sc * (Tsc - Tca) - params.E_el) ...
               / (params.rho_sc * params.A_sc * params.delta_sc * params.cp_sc);
    
    dTca_dt = (h_cd_sc_ca * params.A_ca * (Tsc - Tca) ...
               - h_cd_ca_ct * params.A_ct * (Tca - Tct)) ...
               / (params.rho_ca * params.A_ca * params.delta_ca * params.cp_ca);
    
    dTct_dt = (h_cd_ca_ct * params.A_ct * (Tca - Tct) ...
               - h_cv_ct_wf * params.A_ct * (Tct - Twf)) ...
               / (params.rho_ct * params.A_ct * params.delta_ct * params.cp_ct);
    
    dTwf_dt = (h_cv_ct_wf * params.A_wf * (Tct - Twf) ...
               - params.m_dot_wf * params.cp_wf * (T_w_out - params.T_w_in)) ...
               / (params.rho_wf * params.A_wf * params.delta_wf * params.cp_wf);
    
    dTdt = [dTgs_dt; dTsc_dt; dTca_dt; dTct_dt; dTwf_dt];
end

%% Step 3: Solve for Both Materials
tspan = [0 25200];  % 7 hours
T0 = [30; 30; 30; 30; 30];  % Initial temperatures


fprintf('For Copper Absorber...\n');
[t_cu, T_cu] = ode15s(@(t, T) heatBalanceODE(t, T, params_cu), tspan, T0);

% Solve for Aluminum
fprintf('For Aluminum Absorber...\n');
[t_al, T_al] = ode15s(@(t, T) heatBalanceODE(t, T, params_al), tspan, T0);

%% Step 4: Post-Processing and Comparison
% Convert times to hours
t_hours_cu = t_cu / 3600;
t_hours_al = t_al / 3600;

% Define hour marks (0 to 7 hours)
hours = 0:7;

% Instantaneous temperatures at exact hour marks for both materials
hourly_inst_cu = interp1(t_hours_cu, T_cu, hours, 'linear', 'extrap');
hourly_inst_al = interp1(t_hours_al, T_al, hours, 'linear', 'extrap');

% Layer names
layer_names = {'Tgs', 'Tsc', 'Tca', 'Tct', 'Twf'};

% Print hourly instantaneous values for Copper
fprintf('\nHourly Instantaneous Temperatures (°C) for Copper Absorber:\n');
fprintf('Hour\t');
fprintf('%s\t', layer_names{:});
fprintf('\n');
for i = 1:length(hours)
    fprintf('%d\t', hours(i));
    fprintf('%.2f\t', hourly_inst_cu(i, :));
    fprintf('\n');
end

% Print hourly instantaneous values for Aluminum
fprintf('\nHourly Instantaneous Temperatures (°C) for Aluminum Absorber:\n');
fprintf('Hour\t');
fprintf('%s\t', layer_names{:});
fprintf('\n');
for i = 1:length(hours)
    fprintf('%d\t', hours(i));
    fprintf('%.2f\t', hourly_inst_al(i, :));
    fprintf('\n');
end

% Calculate overall averages for both
avg_cu = mean(T_cu, 1);
avg_al = mean(T_al, 1);

% Print overall averages comparison
fprintf('\nOverall Average Temperatures (°C) Comparison:\n');
fprintf('Layer\tCopper\tAluminum\tDifference (Cu - Al)\n');
for i = 1:5
    diff = avg_cu(i) - avg_al(i);
    fprintf('%s\t%.2f\t%.2f\t%.2f\n', layer_names{i}, avg_cu(i), avg_al(i), diff);
end

% Thermal efficiency for both (using water)
G_interp_cu = interp1(params_cu.time_data, params_cu.G_data, t_cu, 'linear', 'extrap');
eta_th_cu = (params_cu.m_dot_wf * params_cu.cp_wf * (T_cu(:,5) - params_cu.T_w_in)) ./ (G_interp_cu * params_cu.A_gs);

G_interp_al = interp1(params_al.time_data, params_al.G_data, t_al, 'linear', 'extrap');
eta_th_al = (params_al.m_dot_wf * params_al.cp_wf * (T_al(:,5) - params_al.T_w_in)) ./ (G_interp_al * params_al.A_gs);

fprintf('\nAverage Thermal Efficiency Comparison: Cu = %.2f%%, Al = %.2f%%\n', mean(eta_th_cu)*100, mean(eta_th_al)*100);

% Plot Comparison: Temperatures (Combined for Each Layer)
figure;
for i = 1:5
    subplot(3,2,i);  % 5 subplots (adjust as needed)
    plot(t_hours_cu, T_cu(:,i), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Copper');
    hold on;
    plot(t_hours_al, T_al(:,i), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Aluminum');
    xlabel('Time (hours)');
    ylabel('Temperature (°C)');
    title(['Comparison: ' layer_names{i}]);
    legend('Location', 'best');
    grid on;
end

% Plot Efficiency Comparison
figure;
plot(t_hours_cu, eta_th_cu * 100, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Copper');
hold on;
plot(t_hours_al, eta_th_al * 100, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Aluminum');
xlabel('Time (hours)');
ylabel('Thermal Efficiency (%)');
title('Thermal Efficiency Comparison: Cu vs. Al Absorber');
legend('Location', 'best');
grid on;
