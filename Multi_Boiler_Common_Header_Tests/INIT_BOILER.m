% =========================================================================
% INIT_BOILER.M (Improved: energy-consistent steady state & tuned parameters)
% =========================================================================
% - Calibrates k_calibrated from target qdc
% - Computes steady state using BOTH:
%   (i) Q0 = a3*h#*qdc (riser energy, Åström–Bell Eq. 27)
%   (ii) Q0 = qf*(hs - hf) = qs*(hs - hf), with qf=qs at equilibrium (global balance)
% - Loads validation data and pushes all parameters to base workspace
%
% Ref: Åström & Bell (2000), Automatica 36(2000):363–378. 
% =========================================================================

clear; clc;

%% ----------------------- Boiler Initial Conditions -------------------------
p_MPa = 8.7;            % Drum pressure [MPa]
Vwt   = 57.2;           % Water volume [m^3]
a3    = 0.051;          % Steam quality [-]
Vsd   = 4.9;            % Steam volume under drum level [m^3]
qdc   = 1195;           % Natural circulation flow rate [kg/s]
x0    = [Vwt; p_MPa; a3; Vsd];

%% ---------------------------- Geometry & Constants -------------------------
Vr  = 37; Vd = 40; Vdc = 11; Ad = 20; Adc = 0.7; g = 9.81;
b = 0.030; tau_s = 1e-2;
m_metal = 1e6; m_riser = 2e5;

%% ----------------------- Saturated Properties ------------------------------
pbar = p_MPa * 10;                         % [bar]
rho_w = XSteam('rhoL_p', pbar);            % [kg/m^3]
rho_s = XSteam('rhoV_p', pbar);            % [kg/m^3]
h_w   = XSteam('hL_p',  pbar);             % [kJ/kg]
h_s   = XSteam('hV_p',  pbar);             % [kJ/kg]
h_lat = h_s - h_w;                         % latent heat [kJ/kg]

%% ------------------ Riser Steam Volume Fraction ----------------------------
Delta = rho_w - rho_s;
gvoid = a3 * (Delta / max(rho_s,1e-12));
a7    = (rho_w/Delta) * (1 - (rho_s/(rho_w*max(a3,1e-12))) * log(1 + gvoid));

%% ------------------ Loop Friction Coefficient ------------------------------
k_calibrated = (2 * rho_w * Adc * (rho_w - rho_s) * g * a7 * Vr) / max(qdc^2, 1e-12);

%% --------------------------- Steady-State Inputs ---------------------------
Tf0_C = 290;                                % feedwater temperature [C]
hf0   = XSteam('hL_T', Tf0_C);              % feedwater enthalpy [kJ/kg]

Q_kJs = a3 * h_lat * qdc;                   % heat to steam [kJ/s]
Q_MW  = Q_kJs / 1000;                       % [MW]

qs = Q_kJs / max(h_s - hf0, 1e-6);          % steam flow [kg/s]
qf = qs;                                    % feedwater flow [kg/s] (steady-state)

u0 = [Q_MW; qf; qs; hf0];

%% ---------------------------- Load STEAM input ------------------------------
S = readmatrix(fullfile('Digitised Plots', 'Steam Flow Input Step.csv'));
assert(~isempty(S) && size(S,2)>=2, 'Steam Flow Input Step.csv must have [time, q_s] columns.');

tS     = S(:,1);
q_sig  = S(:,2);

tS    = tS(:); q_sig = q_sig(:);
if any(diff(tS)<=0)
    [tS, ix] = unique(tS);
    q_sig = q_sig(ix);
end

% Keep feedwater at steady-state for this single-boiler run
tF     = tS;
qf_sig = qf * ones(size(tF));

steam_input_data.time    = tS;
steam_input_data.signals = struct('values', q_sig, 'dimensions', 1);

%% --------------------- Push variables to base workspace ---------------------
assignin('base','x0',x0);
assignin('base','u0',u0);

assignin('base','b',b);
assignin('base','tau_s',tau_s);
assignin('base','m_metal',m_metal);
assignin('base','m_riser',m_riser);
assignin('base','k_calibrated',k_calibrated);

assignin('base','Vr',Vr);
assignin('base','Vd',Vd);
assignin('base','Vdc',Vdc);
assignin('base','Ad',Ad);
assignin('base','Adc',Adc);

assignin('base','steam_input_data',steam_input_data);

%% ------------------------------- Run Simulink --------------------------------
disp('Boiler Initial Conditions (x0):');
disp(['Vwt = ', num2str(x0(1)), ' m^3']);
disp(['p_MPa = ', num2str(x0(2)), ' MPa']);
disp(['a3 = ', num2str(x0(3)), ' [-]']);
disp(['Vsd = ', num2str(x0(4)), ' m^3']);

disp('Steady-State Inputs (u0):');
disp(['Q_MW = ', num2str(u0(1)), ' MW']);
disp(['qf = ', num2str(u0(2)), ' kg/s']);
disp(['qs = ', num2str(u0(3)), ' kg/s']);
disp(['hf0 = ', num2str(u0(4)), ' kJ/kg']);

simStopTime = num2str(tS(end));
% Update the model name below to your single-boiler Simulink model
sim('Single_Boiler_Common_header.slx', 'StopTime', simStopTime);%