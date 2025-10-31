%% 
%% ================== Two-Boiler Initialization Script (with debug) ==================
clear; clc;

%% ---------------------- User-specified design point ----------------------
p1_MPa = 8.8;   q_s1 = 41;%3  % Boiler 1 drum pressure & flow
p2_MPa = 8.7;   q_s2 = 34.2;%  % Boiler 2 drum pressure & flow

T_f_C  = 290;                  % Feedwater temperature [C]

% Astrom and Bell assumed thermal equilibriumf ro water, steam, and metal
% The approximation is okay as the steam temperature is at 302 deg celsius
% on the simulation

Vwt  = 57.2; 
Vsd0 = 4.9;                    % Nominal steam volume

fprintf('\n=== User Inputs ===\n');
fprintf('Boiler1: p=%.2f MPa, q_s=%.2f kg/s\n', p1_MPa, q_s1);
fprintf('Boiler2: p=%.2f MPa, q_s=%.2f kg/s\n', p2_MPa, q_s2);


%% ---------------------- Geometry & constants -----------------------------
Vr  = 37; Vd = 40; Vdc = 11; Ad = 20; Adc = 0.7;
g = 9.81;
b = 0.030;
tau_s = 1e-2;
m_metal = 1e6;
m_riser = 2e5;

%% ---------------------- Feedwater enthalpy --------------------------------
h_f = XSteam('hL_T', T_f_C);   % kJ/kg
fprintf('\nFeedwater enthalpy h_f = %.2f kJ/kg\n', h_f);

%% ---------------------- Boiler 1 properties -------------------------------
pbar1 = p1_MPa * 10;
rho_w1 = XSteam('rhoL_p', pbar1);
rho_s1 = XSteam('rhoV_p', pbar1);
h_w1   = XSteam('hL_p', pbar1);
h_s1   = XSteam('hV_p', pbar1);
h_lat1 = h_s1 - h_w1;

Q1_kJs = q_s1*h_s1 - q_s1*h_f;  % kJ/s
Q1_MW  = Q1_kJs / 1000;

q_dc = 1195;                     % natural circulation flow [kg/s]
alpha_r1 = fzero(@(a) Q1_kJs - q_dc*a*h_lat1, 0.05);

Vsd1 = Vsd0;

fprintf('\n=== Boiler 1 ===\n');
fprintf('rho_w=%.1f, rho_s=%.1f\n', rho_w1, rho_s1);
fprintf('h_w=%.1f, h_s=%.1f, h_lat=%.1f\n', h_w1, h_s1, h_lat1);
fprintf('Q1=%.2f MW, alpha_r1=%.4f, Vsd1=%.2f m3\n', Q1_MW, alpha_r1, Vsd1);

x0_1 = [Vwt; p1_MPa; alpha_r1; Vsd1];
u0_1 = [Q1_MW; q_s1; q_s1; h_f];

%% ---------------------- Boiler 2 properties -------------------------------
pbar2 = p2_MPa * 10;
rho_w2 = XSteam('rhoL_p', pbar2);
rho_s2 = XSteam('rhoV_p', pbar2);
h_w2   = XSteam('hL_p', pbar2);
h_s2   = XSteam('hV_p', pbar2);
h_lat2 = h_s2 - h_w2;

Q2_kJs = q_s2*h_s2 - q_s2*h_f;   % kJ/s
Q2_MW  = Q2_kJs / 1000;

alpha_r2 = fzero(@(a) Q2_kJs - q_dc*a*h_lat2, 0.05);

Vsd2 = Vsd0;

fprintf('\n=== Boiler 2 ===\n');
fprintf('rho_w=%.1f, rho_s=%.1f\n', rho_w2, rho_s2);
fprintf('h_w=%.1f, h_s=%.1f, h_lat=%.1f\n', h_w2, h_s2, h_lat2);
fprintf('Q2=%.2f MW, alpha_r2=%.4f, Vsd2=%.2f m3\n', Q2_MW, alpha_r2, Vsd2);

x0_2 = [Vwt; p2_MPa; alpha_r2; Vsd2];
u0_2 = [Q2_MW; q_s2; q_s2; h_f];



%% ------------------ Loop Friction Coefficients ----------------------------
Delta1 = rho_w1 - rho_s1;
gvoid1 = alpha_r1 * (Delta1 / rho_s1);
a7_1   = (rho_w1/Delta1) * (1 - (rho_s1/(rho_w1*alpha_r1)) * log(1 + gvoid1));

Delta2 = rho_w2 - rho_s2;
gvoid2 = alpha_r2 * (Delta2 / rho_s2);
a7_2   = (rho_w2/Delta2) * (1 - (rho_s2/(rho_w2*alpha_r2)) * log(1 + gvoid2));

k1 = (2 * rho_w1 * Adc * (rho_w1 - rho_s1) * g * a7_1 * Vr) / max(q_dc^2, 1e-12);
k2 = (2 * rho_w2 * Adc * (rho_w2 - rho_s2) * g * a7_2 * Vr) / max(q_dc^2, 1e-12);

fprintf('\n=== Loop Friction Coefficients ===\n');
fprintf('k1=%.6e, k2=%.6e\n', k1, k2);

%% ---------------------- Push to base workspace ---------------------------
assignin('base','x0_1',x0_1);
assignin('base','x0_2',x0_2);
assignin('base','u0_1',u0_1);
assignin('base','u0_2',u0_2);


assignin('base','b',b);
assignin('base','tau_s',tau_s);
assignin('base','m_metal',m_metal);
assignin('base','m_riser',m_riser);
assignin('base','Vr',Vr);
assignin('base','Vd',Vd);
assignin('base','Vdc',Vdc);
assignin('base','Ad',Ad);
assignin('base','Adc',Adc);
assignin('base','g',g);

assignin('base','k1',k1);
assignin('base','k2',k2);

%% ---------------------- Run Simulink -------------------------------------
simStopTime = '500'; % or based on input signal length
disp('>>> Launching Simulink model...');
sim('Multi_Boiler_Common_header.slx','StopTime',simStopTime);