function [DrumWaterLevel, WaterLevelDeviation] = calculate_level(Vwt, p_MPa, ar, Vsd)
% CALCULATE_LEVEL Computes drum water level and its deviation from nominal.

coder.extrinsic('get_boiler_calcs');

% Geometry and nominal values
Vr   = 37; Vdc = 11; Ad = 20;
Vwt0 = 57.2; Vsd0 = 4.9;
ar0  = 0.051; p0_MPa = 8.66;

% Preallocate doubles
av  = 0.0; av0 = 0.0;

% Call extrinsic and cast
[~, av_tmp]  = get_boiler_calcs(p_MPa, ar);
[~, av0_tmp] = get_boiler_calcs(p0_MPa, ar0);
av  = double(av_tmp);
av0 = double(av0_tmp);

% Water volumes
Vwd0 = Vwt0 - Vdc - (1 - av0) * Vr;
Vwd  = Vwt  - Vdc - (1 - av)  * Vr;

% Levels
l0 = (Vwd0 + Vsd0) / Ad;
DrumWaterLevel = (Vwd + Vsd) / Ad;

% Deviation
WaterLevelDeviation = DrumWaterLevel - l0;
end