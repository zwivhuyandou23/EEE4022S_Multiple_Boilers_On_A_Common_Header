function [mdot_kg_s, h_out_kJkg, dp_valve_MPa, dp_pipe_MPa, ...
          p_after_valve_MPa, p_after_pipe_MPa, mdot_spill_kg_s, mdot_after_spill_kg_s, ...
          regime, Tsat_C, x_out] = ...
    SeriesValvePipe_HEM_DW_SC(p_drum_MPa, x_up, p_hdr_MPa, opening)
% SeriesValvePipe_HEM_DW_SC (IEC-style valve + Darcy–Weisbach pipe + spill)
% Valve uses IEC-style Kv(opening), critical ratio (xt, FL), expansion factor (Fk).
%
% Inputs:
%   p_drum_MPa  upstream drum pressure [MPa]
%   x_up        upstream quality [-] (0..1)
%   p_hdr_MPa   downstream header pressure [MPa]
%   opening     valve opening fraction [0..1]
%
% Outputs:
%   mdot_kg_s             main mass flow through valve+pipe [kg/s]
%   h_out_kJkg            outlet enthalpy [kJ/kg] (HEM, isenthalpic)
%   dp_valve_MPa          pressure drop across valve [MPa]
%   dp_pipe_MPa           pressure drop across pipe [MPa]
%   p_after_valve_MPa     pressure just after valve [MPa]
%   p_after_pipe_MPa      pressure after pipe (header) [MPa]
%   mdot_spill_kg_s       spill mass flow from pipe [kg/s]
%   mdot_after_spill_kg_s main flow minus spill [kg/s]
%   regime                0=no flow, 1=choked, 2=subcritical
%   Tsat_C                upstream saturation temperature [C]
%   x_out                 downstream quality (recalculated at header p)

% ===== Valve constants (IEC-style) =====
Kv_max       = 1.05* 1.8e2;   % [m^3/h] tune to match hardware
Rangeability = 50;             % equal-percentage
Kv_min       = Kv_max / Rangeability;
FL           = 0.90;           % pressure recovery factor
xt           = 0.50;           % choked pressure ratio factor for steam
Fk           = 0.96;           % expansion factor constant for steam

% ===== Pipe & spill (unchanged) =====0.1 * 0.965;
D_m       = 0.1*1.25;       % pipe diameter [m]
L_m       = 50.0;              % pipe length [m]
rough_m   = 4.5e-5;            % absolute roughness [m]
mu_Pa_s   = 1.0e-5;            % dynamic viscosity [Pa·s]
C_SPILL        = 0;%1.0e-2;       % [kg/(s·Pa^0.5)]
p_spill_set_MPa = 8.7;         % spill threshold [MPa]

% ===== Numerics =====
N_iter_max = 25;
alpha      = 0.55;
tol_mdot   = 1e-6;
EPS        = 1e-12;

% ===== Input clamps =====
opening    = max(0.0, min(1.0, double(opening)));
p_drum_MPa = max(0.0, double(p_drum_MPa));
p_hdr_MPa  = max(0.0, double(p_hdr_MPa));
x_up       = max(0.0, min(1.0, double(x_up)));

% ===== Early exit =====
dp_total_Pa = max((p_drum_MPa - p_hdr_MPa) * 1e6, 0.0);
regime = 0.0; Tsat_C = 0.0; x_out = x_up;
if (dp_total_Pa <= 0.0) || (opening <= 0.0)
    [~, h_up_kJkg, Tsat_C] = upstreamProps(p_drum_MPa, x_up);
    mdot_kg_s             = 0.0;
    h_out_kJkg            = h_up_kJkg;
    dp_valve_MPa          = 0.0;
    dp_pipe_MPa           = 0.0;
    p_after_valve_MPa     = p_hdr_MPa;
    p_after_pipe_MPa      = p_hdr_MPa;
    mdot_spill_kg_s       = 0.0;
    mdot_after_spill_kg_s = 0.0;
    x_out                 = recalcQuality(h_out_kJkg, p_hdr_MPa);
    return;
end

% ===== Upstream mixture properties =====
[rho_up, h_up_kJkg, Tsat_C] = upstreamProps(p_drum_MPa, x_up);
h_out_kJkg = h_up_kJkg;  % isenthalpic throttling

% ===== Geometry =====
A_pipe = pi * (D_m/2)^2;

% ===== Iteration state =====
mdot            = max( sqrt(2.0*rho_up*dp_total_Pa) * 1e-3, 0.0 ); % conservative initial guess
dp_pipe_Pa      = 0.0;
dp_valve_Pa     = dp_total_Pa;
p_up_Pa         = p_drum_MPa * 1e6;
p_hdr_Pa        = p_hdr_MPa  * 1e6;
p_after_valve_Pa= p_hdr_Pa;

% ===== Solve series: IEC valve + DW pipe =====
for k = 1:N_iter_max
    % Pipe drop from current mdot
    v_pipe = mdot / max(rho_up * A_pipe, EPS);
    Re     = rho_up * v_pipe * D_m / max(mu_Pa_s, EPS);
    if Re < 2300.0
        f_D = 64.0 / max(Re, 1.0);
    else
        f_D = ( -2.0 * log10( (rough_m/D_m)/3.7 + 5.74 / (max(Re,1.0)^0.9) ) )^(-2);
    end
    f_D = max(f_D, 1e-4);
    dp_pipe_Pa = f_D * (L_m / D_m) * (rho_up * max(v_pipe,0.0)^2 / 2.0);

    % Local downstream pressure right after the valve (before pipe loss)
    p_down_local_Pa = min(max(p_hdr_Pa + dp_pipe_Pa, p_hdr_Pa), p_up_Pa);

    % Valve flow using IEC equations at current p_up/p_dn/opening
    p_up_MPa_loc = p_up_Pa / 1e6;
    p_dn_MPa_loc = p_down_local_Pa / 1e6;

    [mdot_target, dp_eff_bar, regime_loc, Y_loc] = valveIEC_massflow( ...
        p_up_MPa_loc, p_dn_MPa_loc, opening, rho_up, Kv_max, Kv_min, Rangeability, FL, xt, Fk);

    % Update mdot
    mdot_new = alpha * mdot_target + (1.0 - alpha) * mdot;

    % Convergence
    if abs(mdot_new - mdot) <= tol_mdot
        mdot = max(mdot_new, 0.0);
        p_after_valve_Pa = p_down_local_Pa;
        regime = regime_loc;
        break;
    end

    mdot = max(mdot_new, 0.0);
    p_after_valve_Pa = p_down_local_Pa;
    regime = regime_loc;
end

% ===== Spill on pipe =====
if p_drum_MPa > p_spill_set_MPa
    dp_spill_Pa      = (p_drum_MPa - p_spill_set_MPa) * 1e6;
    mdot_spill_kg_s  = C_SPILL * sqrt(max(dp_spill_Pa, 0.0));
else
    mdot_spill_kg_s  = 0.0;
end

% ===== Outputs =====
mdot_kg_s             = max(mdot, 0.0);
mdot_after_spill_kg_s = max(mdot_kg_s - mdot_spill_kg_s, 0.0);

dp_valve_Pa   = max(p_up_Pa - p_after_valve_Pa, 0.0);   % actual remainder drop
dp_valve_MPa  = dp_valve_Pa / 1e6;
dp_pipe_MPa   = dp_pipe_Pa / 1e6;

p_after_valve_MPa = min(max(p_after_valve_Pa / 1e6, p_hdr_MPa), p_drum_MPa);
p_after_pipe_MPa  = p_hdr_MPa;

% Recalculate downstream quality at header pressure
x_out = recalcQuality(h_out_kJkg, p_hdr_MPa);

end % function

% ===== IEC valve mass flow (steam) =====
function [mdot_kg_s, dp_eff_bar, regime, Y] = valveIEC_massflow( ...
    p_up_MPa, p_dn_MPa, opening, rho_up, Kv_max, Kv_min, Rangeability, FL, xt, Fk)

% No-flow check
if (p_dn_MPa >= p_up_MPa) || (opening <= 0.0)
    mdot_kg_s  = 0.0; dp_eff_bar = 0.0; regime = 0.0; Y = 0.0; return;
end

% Upstream/downstream in bar
p_up_bar = 10.0 * max(p_up_MPa, 0.0);
p_dn_bar = 10.0 * max(p_dn_MPa, 0.0);
dp_bar   = max(p_up_bar - p_dn_bar, 0.0);

% Equal-percentage Kv(opening)
Kv_open = min(Kv_max, Kv_min * (Rangeability ^ opening));

% Pressure ratio and critical ratio
x_ratio  = dp_bar / max(p_up_bar, 1e-6);
x_choke  = xt * (FL^2);

if x_ratio > x_choke
    regime = 1.0;                 % choked
    x_eff  = x_choke;             % capacity limited
else
    regime = 2.0;                 % subcritical
    x_eff  = x_ratio;
end

% Effective drop and expansion factor
dp_eff_bar = x_eff * p_up_bar;
Y = max(1.0 - x_eff / (3.0 * Fk), 0.0);

% Volumetric flow (m^3/h), rho_rel uses 1000 kg/m^3 reference
rho_rel = max(rho_up / 1000.0, 1e-6);
Q_m3_h  = Kv_open * Y * sqrt( max(dp_eff_bar, 0.0) / rho_rel );

% Mass flow (kg/s)
mdot_kg_s = (Q_m3_h * rho_up) / 3600.0;
end

% ===== Upstream props via XSteam =====
function [rho_up, h_up_kJkg, Tsat_C] = upstreamProps(p_up_MPa, x_up)
p_bar = 10.0 * p_up_MPa;
Tsat_C = xs_p('Tsat_p', p_bar);
hf = xs_p('hL_p', p_bar); hg = xs_p('hV_p', p_bar);
vf = xs_p('vL_p', p_bar); vg = xs_p('vV_p', p_bar);
v_mix = max((1.0 - x_up)*vf + x_up*vg, 1e-9);
rho_up = 1.0 / v_mix;
h_up_kJkg = (1.0 - x_up)*hf + x_up*hg;
end

% ===== Downstream quality from enthalpy at header p =====
function x_out = recalcQuality(h_out_kJkg, p_hdr_MPa)
p_out_bar = 10.0 * max(p_hdr_MPa, 0.0);
hf_out = xs_p('hL_p', p_out_bar); hg_out = xs_p('hV_p', p_out_bar);
if h_out_kJkg <= hf_out
    x_out = 0.0;
elseif h_out_kJkg >= hg_out
    x_out = 1.0;
else
    x_out = (h_out_kJkg - hf_out) / max(hg_out - hf_out, 1e-9);
end
end

% ===== XSteam typed wrapper =====
function val = xs_p(token, pbar)
coder.extrinsic('XSteam');
val = 0.0;
tmp = feval('XSteam', token, pbar);
if ~isempty(tmp), val = double(tmp); end
if ~isfinite(val), val = 0.0; end
end