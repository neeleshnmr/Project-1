clear; clc; close all;
params = struct();


params.A_gs = 0.5562;      
params.A_sc = 0.4903;      
params.A_Ala = 0.4903;      
params.A_ct = 0.4903 * 0.1; 
params.A_wf = 0.4903 * 0.1; 


params.delta_gs = 0.003;   
params.delta_sc = 0.0003;  
params.delta_Ala = 0.0015;  
params.delta_ct = 0.01;    
params.delta_wf = 0.01;    

params.rho_gs = 2500;      
params.rho_sc = 2330;      
params.rho_Ala = 2700;      
params.rho_ct = 8960;      
params.rho_wf = 1000;      


params.cp_gs = 750;        
params.cp_sc = 700;        
params.cp_Ala = 900;        
params.cp_ct = 385;        
params.cp_wf = 4180;       

% Flow rates and inlets
params.m_dot_wf = 0.05;    
params.T_w_in = 30;        


params.time_data = [0, 3600, 7200, 10800, 14400, 18000, 21600, 25200];  % Seconds
params.G_data = [200, 600, 900, 1100, 1000, 800, 500, 300];  % W/m²
params.T_am_data = [25, 28, 30, 32, 31, 30, 29, 28];  % °C


params.alpha_gs = 0.15;    
params.alpha_sc = 0.9;    
params.E_el = 30;          

% Parameters for dynamic h calculations (adapted, air-related removed)
params.epsilon_g = 0.9;    % Glass emissivity
params.sigma = 5.67e-8;    % Stefan-Boltzmann constant [W/m²·K⁴]
params.Vw = 1.3;             % Wind speed [m/s] for conv coeff
params.K_g = 1.8;         
params.K_sc = 148;         
params.K_Ala = 237;         
params.K_ct = 385;         




function dTdt = heatBalanceODE(t, T, params)
    
    Tgs = T(1); Tsc = T(2); TAla = T(3); Tct = T(4); Twf = T(5);
    
    % Interpolate time-varying inputs
    G = interp1(params.time_data, params.G_data, t, 'linear', 'extrap');
    T_am = interp1(params.time_data, params.T_am_data, t, 'linear', 'extrap');
    T_s = 0.0552 * T_am^1.5;  % Sky temp
    
   
    % Radiative: glass to sky (Eq. 13)
    h_rd_gs_s = params.epsilon_g * params.sigma * (Tgs^2 + T_s^2) * (Tgs + T_s);
    
    % Convective: glass to ambient (Eq. 15)
    h_cv_gs_am = 2.8 + 3 * params.Vw;
    
    % Conductive: between layers (Eq. 16, etc.)
    h_cd_gs_sc = 1 / (params.delta_gs / params.K_g + params.delta_sc / params.K_sc);
    h_cd_sc_Ala = 1 / (params.delta_sc / params.K_sc + params.delta_Ala / params.K_Ala);
    h_cd_ca_ct = 1 / (params.delta_Ala / params.K_Ala + params.delta_ct / params.K_ct);
    
   
    h_cv_ct_wf = 2745;  %can nusselt no be used?
    
  
    Egs = params.alpha_gs * G * params.A_gs;  % Eq. 6
    Esc = params.alpha_sc * G * params.A_sc;  % Eq. 7
    
    
    T_w_out = Twf;
    
    
    % Equation (6): dTgs/dt
    dTgs_dt = (Egs - h_rd_gs_s * params.A_gs * (Tgs - T_s) ...
               - h_cv_gs_am * params.A_gs * (Tgs - T_am) ...
               - h_cd_gs_sc * params.A_gs * (Tgs - Tsc)) ...
               / (params.rho_gs * params.A_gs * params.delta_gs * params.cp_gs);
    
    % Equation (7): dTsc/dt
    dTsc_dt = (Esc + h_cd_gs_sc * params.A_sc * (Tgs - Tsc) ...
               - h_cd_sc_Ala * params.A_sc * (Tsc - TAla) - params.E_el) ...
               / (params.rho_sc * params.A_sc * params.delta_sc * params.cp_sc);
    
    % Adapted Equation (8): dTca/dt 
    dTca_dt = (h_cd_sc_Ala * params.A_Ala * (Tsc - TAla) ...
               - h_cd_ca_ct * params.A_ct * (TAla - Tct)) ...
               / (params.rho_Ala * params.A_Ala * params.delta_Ala * params.cp_Ala);
    
    % Adapted Equation (11): dTct/dt 
    dTct_dt = (h_cd_ca_ct * params.A_ct * (TAla - Tct) ...
               - h_cv_ct_wf * params.A_ct * (Tct - Twf)) ...
               / (params.rho_ct * params.A_ct * params.delta_ct * params.cp_ct);
    
    % Equation (12): dTwf/dt
    dTwf_dt = (h_cv_ct_wf * params.A_wf * (Tct - Twf) ...
               - params.m_dot_wf * params.cp_wf * (T_w_out - params.T_w_in)) ...
               / (params.rho_wf * params.A_wf * params.delta_wf * params.cp_wf);
    
    % Return derivatives (now 5)
    dTdt = [dTgs_dt; dTsc_dt; dTca_dt; dTct_dt; dTwf_dt];
end


tspan = [0 25200];  % 7 hours (25200 seconds)
T0 = [30; 30; 30; 30; 30];  % Initial temperatures = ambient [°C] (5 layers)

[t, T] = ode15s(@(t, T) heatBalanceODE(t, T, params), tspan, T0);


% Useful thermal power (water, adapted Eq. 2)
% Qu_wf = params.m_dot_wf * params.cp_wf * (T(:,5) - params.T_w_in);  % [W]

% Thermal efficiency (Eq. 1) - With interpolated G for dynamic efficiency
% G_interp = interp1(params.time_data, params.G_data, t, 'linear', 'extrap');
% eta_th = Qu_wf ./ (G_interp * params.A_gs);  % Fraction

% Convert time to hours
t_hours = t / 3600;

% Define hour marks (0 to 7 hours)
hours = 0:7;

% Instantaneous temperatures at exact hour marks (using interpolation)
hourly_inst = interp1(t_hours, T, hours, 'linear', 'extrap');  % Rows: hours, Columns: layers

% Layer names for printing
layer_names = {'Tgs', 'Tsc', 'TAla', 'Tct', 'Twf'};

% Print instantaneous hourly values for all parts
fprintf('\nHourly Instantaneous Temperatures (°C) for All Parts:\n');
fprintf('Hour\t');
fprintf('%s\t', layer_names{:});
fprintf('\n');
for i = 1:length(hours)
    fprintf('%d\t', hours(i));
    fprintf('%.2f\t', hourly_inst(i, :));
    fprintf('\n');
end

% Calculate overall averages
avg_Tgs = mean(T(:,1));   % Glass surface
avg_Tsc = mean(T(:,2));   % Solar cells
avg_TAla = mean(T(:,3));   % Al absorber
avg_Tct = mean(T(:,4));   % Copper tube
avg_Twf = mean(T(:,5));   % Water flow

% Print overall averages
fprintf('\nOverall Average Temperatures (°C):\n');
fprintf('- Glass Surface (Tgs): %.2f\n', avg_Tgs);
fprintf('- Solar Cells (Tsc): %.2f\n', avg_Tsc);
fprintf('- Aluminium Absorber (TAla): %.2f\n', avg_TAla);
fprintf('- Copper Tube (Tct): %.2f\n', avg_Tct);
fprintf('- Water Flow (Twf): %.2f\n', avg_Twf);
%fprintf('\nOverall Average Thermal Efficiency (Water): %.2f%%\n', mean(eta_th) * 100);

% Combined time-varying graph for all temperatures
figure(1);
plot(t_hours, T, 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend(layer_names, 'Location', 'best');
title('Time-Varying Temperature Profiles (Water-Only System)');
grid on;

% Thermal efficiency graph
% figure(2);
% plot(t_hours, eta_th * 100, 'LineWidth', 1.5, 'Color', 'r');  % Convert to %
% xlabel('Time (hours)');
% ylabel('Thermal Efficiency (%)');
% title('Time-Varying Water Thermal Efficiency');
% grid on;
