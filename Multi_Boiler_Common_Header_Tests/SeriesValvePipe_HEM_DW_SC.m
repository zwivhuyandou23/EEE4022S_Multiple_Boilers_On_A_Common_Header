function [mdot_kg_s, h_out_kJkg, dp_valve_MPa, dp_pipe_MPa, ...
          p_after_valve_MPa, p_after_pipe_MPa, mdot_spill] = ...
    SeriesValvePipe_HEM_DW_SC(p_drum_MPa, x_up, p_hdr_MPa, opening)
% SeriesValvePipe_HEM_DW_SC (Kv-based, with pipe spill)
% Computes flow through a valve in series with a pipe. 
% Adds a spill term on the pipe (e.g. leakage/relief) proportional to Δp.
%
% Inputs:
%   p_drum_MPa   upstream (boiler drum) pressure [MPa]
%   x_up         upstream quality [-] (0..1)
%   p_hdr_MPa    downstream (header) pressure [MPa]
%   opening      valve opening fraction [0..1]
%
% Outputs:
%   mdot_kg_s        main mass flow through valve+pipe [kg/s]
%   h_out_kJkg       outlet enthalpy [kJ/kg] (HEM, isenthalpic)
%   dp_valve_MPa     pressure drop across valve [MPa]
%   dp_pipe_MPa      pressure drop across pipe [MPa]
%   p_after_valve_MPa pressure just after valve [MPa]
%   p_after_pipe_MPa  pressure after pipe (header) [MPa]
%   mdot_spill       spill mass flow from pipe [kg/s]

% ===== Valve constants (Kv-based) =====
CD_OUT  = 0.70;       % discharge coefficient [-]
KV_MAX  = 90.0e-2;       % max valve coefficient [m^3/h] @ water 20°C
R_EQP   = 50.0;       % equal-percentage rangeability [-]

% ===== Pipe constants =====
D_m      = 0.1*0.975;      % pipe diameter [m]
L_m      = 50.0;      % pipe length [m]
rough_m  = 4.5e-5;    % absolute roughness [m]
mu_Pa_s  = 1.0e-5;    % dynamic viscosity [Pa·s]

% ===== Spill constants (pipe only) =====
C_SPILL  = 0;%1.0e-3;    % spill coefficient [kg/(s·Pa^0.5)]
p_spill_set = 8.0;    % MPa, threshold for spill onset

% ===== Numerics =====
N_iter   = 8;
alpha    = 0.5;

% ===== Safety clamps =====
opening     = max(0.0, min(1.0, double(opening)));
p_drum_MPa  = double(p_drum_MPa);
p_hdr_MPa   = double(p_hdr_MPa);
x_up        = max(0.0, min(1.0, double(x_up)));

% ===== Total ΔP available =====
dp_total_Pa = max((p_drum_MPa - p_hdr_MPa) * 1e6, 0.0);

if (dp_total_Pa <= 0.0) || (opening <= 0.0)
    mdot_kg_s        = 0.0;
    h_out_kJkg       = 0.0;
    dp_valve_MPa     = 0.0;
    dp_pipe_MPa      = 0.0;
    p_after_valve_MPa= p_hdr_MPa;
    p_after_pipe_MPa = p_hdr_MPa;
    mdot_spill       = 0.0;
    return;
end

% ===== Upstream mixture properties at drum pressure =====
p_bar = 10.0 * p_drum_MPa;
hf = xs('hL_p', p_bar); hg = xs('hV_p', p_bar);
vf = xs('vL_p', p_bar); vg = xs('vV_p', p_bar);
v_up  = max((1.0 - x_up)*vf + x_up*vg, 1.0e-9);
rho_up = 1.0 / v_up;
h_up_kJkg = (1.0 - x_up)*hf + x_up*hg;

% Assume isenthalpic HEM
h_out_kJkg = h_up_kJkg;

% ===== Effective Kv (equal-percentage trim) =====
if R_EQP > 1.0
    f = (R_EQP.^opening - 1.0) / (R_EQP - 1.0);
else
    f = opening;
end
Kv_eff = max(0.0, KV_MAX * f);
Aeff   = Kv_eff / (161.0 * CD_OUT);

% ===== Pipe geometry =====
A_pipe = pi * (D_m/2)^2;

% ===== Iterative series solve for mdot =====
mdot = 41;%Aeff * sqrt(2.0 * rho_up * dp_total_Pa) * 0.5;
dp_pipe_Pa  = 0.0;
dp_valve_Pa = dp_total_Pa;
p_after_valve_Pa = p_hdr_MPa*1e6;

for k = 1:N_iter
    % Pipe drop (Darcy–Weisbach)
    v_pipe = mdot / (rho_up * A_pipe);
    Re     = rho_up * v_pipe * D_m / mu_Pa_s;
    if Re < 2300.0
        f_D = 64.0 / max(Re, 1.0e-6);
    else
        f_D = ( -2.0 * log10( (rough_m/D_m)/3.7 + 5.74 / (max(Re,1.0e-6)^0.9) ) )^(-2);
    end
    dp_pipe_Pa = f_D * (L_m / D_m) * (rho_up * v_pipe^2 / 2.0);

    % Downstream local pressure after pipe
    p_down_local_Pa = p_hdr_MPa * 1e6 + dp_pipe_Pa;

    % Valve drop is the remainder
    dp_valve_Pa = max(p_drum_MPa*1e6 - p_down_local_Pa, 0.0);

    % Update mdot from valve (HEM orifice)
    mdot_target = CD_OUT * Aeff * sqrt(max(2.0 * rho_up * dp_valve_Pa, 0.0));

    % Relax for stability
    mdot = alpha * mdot_target + (1.0 - alpha) * mdot;

    % Capture pressure after valve
    p_after_valve_Pa = p_down_local_Pa;
end

% ===== Spill flow (pipe aspect only) =====
if p_drum_MPa > p_spill_set
    dp_spill_Pa = (p_drum_MPa - p_spill_set) * 1e6;
    mdot_spill  = C_SPILL * sqrt(dp_spill_Pa);
else
    mdot_spill  = 0.0;
end

% ===== Final outputs =====
mdot_kg_s        = mdot;
dp_valve_MPa     = dp_valve_Pa / 1e6;
dp_pipe_MPa      = dp_pipe_Pa  / 1e6;
p_after_valve_MPa= p_after_valve_Pa / 1e6;
p_after_pipe_MPa = p_hdr_MPa;
end

% ===== Local helper =====
function val = xs(token, pbar)
coder.extrinsic('XSteam');
val = 0.0;
tmp = feval('XSteam', token, pbar);
if ~isempty(tmp), val = double(tmp); end
if ~isfinite(val), val = 0.0; end
end