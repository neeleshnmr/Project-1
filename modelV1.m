clear; clc; close all;

params = struct();

params.A_gs = 0.5562;      
params.A_sc = 0.4903;      
params.A_ca = 0.4903;      
params.A_af = 0.4903;      
params.A_cfi = 0.4903;     
params.A_ct = 0.4903 * 0.1; 
params.A_wf = 0.4903 * 0.1; 


params.delta_gs = 0.003;   
params.delta_sc = 0.0003;  
params.delta_ca = 0.0015;  
params.delta_af = 0.05;    
params.delta_cfi = 0.0015; 
params.delta_ct = 0.01;    
params.delta_wf = 0.01;    


params.rho_gs = 2500;      
params.rho_sc = 2330;      
params.rho_ca = 8960;      
params.rho_af = 1.2;       
params.rho_cfi = 8960;     
params.rho_ct = 8960;      
params.rho_wf = 1000;      


params.cp_gs = 750;        
params.cp_sc = 700;        
params.cp_ca = 385;        
params.cp_af = 1005;       
params.cp_cfi = 385;       
params.cp_ct = 385;        
params.cp_wf = 4180;       


params.h_rd_gs_s = 5;      
params.h_cv_gs_am = 10;    
params.h_cd_gs_sc = 200;   
params.h_cd_sc_ca = 200;  
params.h_cd_ca_ct = 300;   
params.h_cd_ca_cfi = 300;  
params.h_cv_ca_af = 20;    
params.h_cv_cfi_af = 20;   
params.h_cv_ct_af = 20;    
params.h_cv_ct_wf = 100;   


params.m_dot_af = 0.05;    % air flow rate
params.m_dot_wf = 0.02;    % water flow rate
params.T_a_in = 30;        % inlet air temp
params.T_w_in = 30;        % inlet water telp

params.time_data = [0, 3600, 7200, 10800, 14400, 18000, 21600, 25200];
params.G_data = [200, 600, 900, 1100, 1000, 800, 500, 300];
params.T_am_data = [25, 28, 30, 32, 31, 30, 29, 28];

params.G = 900;            
params.alpha_gs = 0.05;    
params.alpha_sc = 0.9;     
params.T_am = 30;          
params.T_s = 0.0552 * params.T_am^1.5;  


params.E_el = 30;          %electrical output, assumed


function dTdt = heatBalanceODE(t, T, params)
   
    Tgs = T(1); Tsc = T(2); Tca = T(3); Taf = T(4); Tcfi = T(5); Tct = T(6); Twf = T(7);
    G = interp1(params.time_data, params.G_data, t, 'linear', 'extrap');
    T_am = interp1(params.time_data, params.T_am_data, t, 'linear', 'extrap');
    T_s = 0.0552 * T_am^1.5;  
    
    
    Egs = params.alpha_gs * params.G * params.A_gs;  % using assumed 900 w/m2 val
    Esc = params.alpha_sc * params.G * params.A_sc;  % same as above reeason
    
    % Assume outlet temps approx. equal to fluid temps for advection (simplified)
    T_a_out = Taf;  % Air outlet
    T_w_out = Twf;  % Water outlet
    
    % Equation (6): dTgs/dt
    dTgs_dt = (Egs - params.h_rd_gs_s * params.A_gs * (Tgs - T_s) ...
               - params.h_cv_gs_am * params.A_gs * (Tgs - params.T_am) ...
               - params.h_cd_gs_sc * params.A_gs * (Tgs - Tsc)) ...
               / (params.rho_gs * params.A_gs * params.delta_gs * params.cp_gs);
    
    % Equation (7): dTsc/dt
    dTsc_dt = (Esc + params.h_cd_gs_sc * params.A_sc * (Tgs - Tsc) ...
               - params.h_cd_sc_ca * params.A_sc * (Tsc - Tca) - params.E_el) ...
               / (params.rho_sc * params.A_sc * params.delta_sc * params.cp_sc);
    
    % Equation (8): dTca/dt
    dTca_dt = (params.h_cd_sc_ca * params.A_ca * (Tsc - Tca) ...
               - params.h_cd_ca_ct * params.A_ct * (Tca - Tct) ...
               - params.h_cd_ca_cfi * params.A_cfi * (Tca - Tcfi) ...
               - params.h_cv_ca_af * params.A_ca * (Tca - Taf)) ...
               / (params.rho_ca * params.A_ca * params.delta_ca * params.cp_ca);
    
    % Equation (9): dTaf/dt
    dTaf_dt = (params.h_cv_ca_af * params.A_ca * (Tca - Taf) ...
               + params.h_cv_cfi_af * params.A_cfi * (Tcfi - Taf) ...
               + params.h_cv_ct_af * params.A_ct * (Tct - Taf) ...
               - params.m_dot_af * params.cp_af * (T_a_out - params.T_a_in)) ...
               / (params.rho_af * params.A_af * params.delta_af * params.cp_af);
    
    % Equation (10): dTcfi/dt
    dTcfi_dt = (params.h_cd_ca_cfi * params.A_cfi * (Tca - Tcfi) ...
                - params.h_cv_cfi_af * params.A_cfi * (Tcfi - Taf)) ...
                / (params.rho_cfi * params.A_cfi * params.delta_cfi * params.cp_cfi);
    
    % Equation (11): dTct/dt
    dTct_dt = (params.h_cd_ca_ct * params.A_ct * (Tca - Tct) ...
               - params.h_cv_ct_af * params.A_ct * (Tct - Taf) ...
               - params.h_cv_ct_wf * params.A_ct * (Tct - Twf)) ...
               / (params.rho_ct * params.A_ct * params.delta_ct * params.cp_ct);
    
    % Equation (12): dTwf/dt
    dTwf_dt = (params.h_cv_ct_wf * params.A_wf * (Tct - Twf) ...
               - params.m_dot_wf * params.cp_wf * (T_w_out - params.T_w_in)) ...
               / (params.rho_wf * params.A_wf * params.delta_wf * params.cp_wf);
    
    % Return derivatives
    dTdt = [dTgs_dt; dTsc_dt; dTca_dt; dTaf_dt; dTcfi_dt; dTct_dt; dTwf_dt];
end


tspan = [0 25200];  
T0 = [30; 30; 30; 30; 30; 30; 30]; 

[t, T] = ode45(@(t, T) heatBalanceODE(t, T, params), tspan, T0);


Qu_af = params.m_dot_af * params.cp_af * (T(:,4) - params.T_a_in);  % [W]

% Thermal efficiency (Eq. 1)
eta_th = Qu_af ./ (params.G * params.A_gs);

t_hours = t / 3600;
hours = 0:7;
hourly_inst = interp1(t_hours, T, hours, 'linear', 'extrap'); 
layer_names = {'Tgs', 'Tsc', 'Tca', 'Taf', 'Tcfi', 'Tct', 'Twf'};


fprintf('\nInstantaneous Temperatures (°C) at Each Hour:\n');
fprintf('Hour\t');
fprintf('%s\t', layer_names{:});
fprintf('\n');
for i = 1:length(hours)
    fprintf('%d\t', hours(i));
    fprintf('%.2f\t', hourly_inst(i, :));
    fprintf('\n');
end

% Plot temperatures
% figure;
% plot(t/3600, T);
% xlabel('Time (hours)');
% ylabel('Temperature (°C)');
% legend('Tgs', 'Tsc', 'Tca', 'Taf', 'Tcfi', 'Tct', 'Twf');
% title('Temperature Profiles in PV/T System');

% Plot thermal efficiency
% figure;
% plot(t/3600, eta_th * 100);  % Convert to %
% xlabel('Time (hours)');
% ylabel('Thermal Efficiency (%)');
% title('Air Thermal Efficiency');

% figure(1);
% plot(t_hours, T, 'LineWidth', 1.5);
% xlabel('Time (hours)');
% ylabel('Temperature (°C)');
% legend(layer_names, 'Location', 'best');
% title('Time-Varying Temperature Profiles for All Layers');
% grid on;


% Display average values

avg_Tgs = mean(T(:,1));   
avg_Tsc = mean(T(:,2));   
avg_Tca = mean(T(:,3));   
avg_Taf = mean(T(:,4));   
avg_Tcfi = mean(T(:,5));  
avg_Tct = mean(T(:,6));   
avg_Twf = mean(T(:,7));
fprintf('Average Temperatures (°C):\n');
fprintf('- Glass Surface (Tgs): %.2f\n', avg_Tgs);
fprintf('- Solar Cells (Tsc): %.2f\n', avg_Tsc);
fprintf('- Copper Absorber (Tca): %.2f\n', avg_Tca);
fprintf('- Airflow (Taf): %.2f\n', avg_Taf);
fprintf('- Copper Fins (Tcfi): %.2f\n', avg_Tcfi);
fprintf('- Copper Tube (Tct): %.2f\n', avg_Tct);
fprintf('- Water Flow (Twf): %.2f\n', avg_Twf);
% fprintf('Average Thermal Efficiency (Air): %.2f%%\n', mean(eta_th) * 100);
