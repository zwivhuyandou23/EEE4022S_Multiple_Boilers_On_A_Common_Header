function q_f = feedwater_valve_boiler2(m)
% Feedwater control valve model (MPa input), using proper Kvs liquid formula
% K_vs: m^3/hr per sqrt(bar), water at ~20 C
% m:    normalized opening [0..1]
% DeltaP_MPa: pressure drop [MPa]
% rho:  density [kg/m^3]

% Design parameters (tune as needed)
K_vs        = 144.7547;    % set ~100 for ~40 kg/s at m~0.75, DeltaP=0.5 MPa, T~290 C
T_f_C       = 290;    % feedwater temperature [C]
DeltaP_MPa  = 0.4;    % pressure drop [MPa]

% Density (better: use T & P; here saturated liquid at T)
rho = XSteam('rhoL_T', T_f_C);   % kg/m^3

% Convert MPa -> bar
DeltaP_bar = DeltaP_MPa * 10;

% Relative density SG
SG = rho / 1000;

% Volumetric flow (m^3/hr) using Kvs liquid formula
Q_m3_hr = K_vs * max(0,min(1,m)) * sqrt( max(DeltaP_bar,0) / max(SG,eps) );

% Mass flow (kg/s)
q_f = (rho * Q_m3_hr) / 3600; % 34.2 target at steady state
end