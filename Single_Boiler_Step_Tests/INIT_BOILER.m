% =========================================================================
% INIT_BOILER.M — Calibrated Initialization for Åström & Bell Drum-Boiler
% (No plotting; prepares From Workspace timeseries inputs for Simulink)
% =========================================================================
clear; clc;

fprintf('====================================================\n');
fprintf('    INITIALIZING & CALIBRATING BOILER STEADY-STATE\n');
fprintf('====================================================\n\n');

%% ------------------ PART 1: Define Target State ------------------
p0_MPa      = 8.5;        % Drum pressure [MPa]
Vwt0        = 57.2;       % Water volume [m^3]
a30         = 0.051;      % Riser outlet steam quality [-]
Vsd0        = 4.9;        % Steam volume under drum level [m^3]
qdc0_target = 1195;       % Circulation flow [kg/s]

x0 = [Vwt0; p0_MPa; a30; Vsd0];

%% ------------------ PART 2: Physical Parameters ------------------
Vr       = 37; Vd = 40; Vdc = 11; Ad = 20; Adc = 0.7;
g        = 9.81;
b      = 0.030e2;    % Empirical coefficient in q_4d (dimensionless)
tau_s  = 1e-2;    % Steam residence time in drum [s]

% Metal masses (dominant in pressure dynamics)
m_metal = 0.3*(1e6);  % Total boiler metal mass [kg]
m_riser = 2.0e5;  % Riser metal mass [kg]

%% ------------------ PART 3: Steam Properties ------------------
p0_bar = p0_MPa * 10;
props = get_steam_props(p0_bar);

rho_w   = props.rhoL; rho_s = props.rhoV;
h_w     = props.hL;   h_s   = props.hV;
v_w     = props.vL;   v_s   = props.vV;

h_fg    = h_s - h_w;
av0     = avbar_corrected(a30, v_w, v_s);

%% ------------------ PART 4: Calibrate Friction Factor ------------------
circ_numerator = 2 * rho_w * Adc^2 * (rho_w - rho_s) * g * av0 * Vr;
k1   = circ_numerator / (qdc0_target^2);

%% ------------------ PART 5: Compute Consistent Inputs ------------------
Tf0_C = 290;    %                      
hf0   = XSteam('hL_T', Tf0_C);       

qf0   = a30 * qdc0_target;           
qs0   = qf0;                         
Q0_MW = qf0 * (h_s - hf0) / 1000;    

u0 = [Q0_MW; qf0; qs0; hf0];

%% ------------------ PART 6: Load Digitised Reference Data ------------------
fprintf('\n--- Loading Digitised Reference Data ---\n');

dataFolder = fullfile(pwd,'Digitised Plots');

files = { ...
    'Heat_Step_Drum_Pressure.csv', ...
    'Heat_Step_Drum_water_Level.csv', ...
    'Heat_Step_Steam_Quality.csv', ...
    'Heat_Step_Total_Water_Volume.csv', ...
    'Heat_Step_Volume_of_Steam_In_Drum.csv', ...
    'Steam Flow Input Step.csv', ...
    'Steam_Flow_Step_Drum_pressure.csv', ...
    'Steam_Flow_Step_Drum_water_Level.csv', ...
    'Steam_Flow_Step_Steam_Quality.csv', ...
    'Steam_Flow_Step_Total_Water_Volume.csv', ...
    'Steam pressure.csv', ...
    'Water level deviation.csv'};

digitised = struct();

for i = 1:numel(files)
    fpath = fullfile(dataFolder, files{i});
    if isfile(fpath)
        raw = readmatrix(fpath);
        t = raw(:,1); y = raw(:,2);
        [~,fname,~] = fileparts(files{i});
        field = matlab.lang.makeValidName(fname);
        digitised.(field).t = t;
        digitised.(field).y = y;
        % Also create timeseries for Simulink
        ts = timeseries(y,t);
        assignin('base',[field '_ts'],ts);
        fprintf('Loaded %-40s (%d points)\n', files{i}, length(t));
    else
        warning('File not found: %s', files{i});
    end
end

assignin('base','digitised',digitised);

%% ------------------ PART 7: Push Calibration Variables ------------------
assignin('base','x0',x0);
assignin('base','u0',u0);
assignin('base','b',b);

assignin('base','k', k1);

assignin('base','Vr',Vr);
assignin('base','Vd',Vd);
assignin('base','Vdc',Vdc);
assignin('base','Ad',Ad);
assignin('base','Adc',Adc);
assignin('base','g',g);
assignin('base','m_metal',m_metal);
assignin('base','m_riser',m_riser);
assignin('base','tau_s',tau_s);
assignin('base','hf0',hf0);

fprintf('\n--- All calibrated variables and timeseries loaded into workspace ---\n');

%% ------------------ HELPER FUNCTIONS ------------------
function props = get_steam_props(p_bar)
    props.rhoL = XSteam('rhoL_p', p_bar);
    props.rhoV = XSteam('rhoV_p', p_bar);
    props.hL   = XSteam('hL_p',  p_bar);
    props.hV   = XSteam('hV_p',  p_bar);
    props.vL   = XSteam('vL_p',  p_bar);
    props.vV   = XSteam('vV_p',  p_bar);
end

function av = avbar_corrected(a3, v_w, v_s)
    denom = v_s*(1 - a3) + v_w*a3;
    term1 = v_w / max(denom, 1e-12);
    gamma = a3 * (v_s - v_w) / max(denom, 1e-12);
    if abs(gamma) < 1e-9
        log_term = 1 - gamma/2 + gamma^2/3;
    else
        log_term = log(1 + gamma) / gamma;
    end
    av = term1 * log_term;
    av = min(max(av,0.0),1.0);
end

disp('Initial States (x0): [Vwt0, p0_MPa, a30, Vsd0]');
disp(x0);

disp('Initial Inputs (u0): [Q0_MW, qf0, qs0, hf0]');
disp(u0);
%sim("Single_Boiler_Model.slx")