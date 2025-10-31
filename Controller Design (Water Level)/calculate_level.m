function [DrumWaterLevel, WaterLevelDeviation, av_out, WaterLevelContribution, SteamLevelContribution] = calculate_level(Vwt, p_MPa, ar, Vsd)
% CALCULATE_LEVEL Computes drum water level and its deviation from nominal.
% 
% Inputs:
%   Vwt   - Total water volume in the system [m^3]
%   p_MPa - Drum pressure [MPa]
%   ar    - Steam quality at riser outlet [-]
%   Vsd   - Steam volume under drum level [m^3]
%
% Outputs:
%   DrumWaterLevel        - Actual drum water level [m]
%   WaterLevelDeviation   - Deviation from nominal level [m]
%   av_out                - Average riser void fraction [-]
%   WaterLevelContribution- Contribution to level from water volume [m]



coder.extrinsic('get_boiler_calcs');

% --- Fixed geometry and nominal operating point values ---
Vr = 37;       % Riser volume [m^3]
Vdc = 11;      % Downcomer volume [m^3]
Ad = 20;       % Drum cross-sectional area [m^2]
Vwt0 = 57.16;  % Nominal total water volume [m^3]
Vsd0 = 4.9;    % Nominal steam volume under drum level [m^3]
ar0 = 0.051;   % Nominal riser outlet steam quality [-]
p0_MPa = 8.5;  % Nominal drum pressure [MPa]

p0_MPa = 8.66;          % [MPa]
Vwt0   = 57.16;          % [m^3]
a30    = 0.051;         % [-]
Vsd0   = 4.9;           % [m^3]
qdc0_target = 1195;     % [kg/s] primary calibration target


DrumWaterLevel = 0.0;
WaterLevelDeviation = 0.0;
av_out = 0.0;
WaterLevelContribution = 0.0;
SteamLevelContribution = 0.0;


props = struct('rhoL',0.0, 'rhoV',0.0, 'hL',0.0, 'hV',0.0, 'vL',0.0, 'vV',0.0, ...
               'Tsat',0.0, ...
               'drhoL_dpbar',0.0, 'drhoV_dpbar',0.0, 'dhL_dpbar',0.0, 'dhV_dpbar',0.0, ...
               'dvL_dpbar',0.0, 'dvV_dpbar',0.0, 'dTs_dpbar',0.0);
props0 = props;  % Copy structure for nominal state


av = 0.0;
av0 = 0.0;


[props, av] = get_boiler_calcs(p_MPa, ar);
v_w = props.vL;  % Specific volume of water
v_s = props.vV;  % Specific volume of steam


[props0, av0] = get_boiler_calcs(p0_MPa, ar0);
v_w0 = props0.vL;
v_s0 = props0.vV;


Vwd0 = Vwt0 - Vdc - (1 - av0) * Vr;  % Nominal drum water volume
Vwd  = Vwt  - Vdc - (1 - av)  * Vr;  % Current drum water volume


l0 = (Vwd0 + Vsd0) / Ad;             % Nominal drum level
DrumWaterLevel = (Vwd + Vsd) / Ad;   % Current drum level


WaterLevelDeviation =  (DrumWaterLevel-l0) ;  % Negative sign: deviation below nominal is positive

% --- Compute contributions to level from water and steam volumes ---
WaterLevelContribution = (Vwd / Ad) - (Vwd0 / Ad);
SteamLevelContribution = (Vsd / Ad) - (Vsd0 / Ad);


av_out = av;
end

