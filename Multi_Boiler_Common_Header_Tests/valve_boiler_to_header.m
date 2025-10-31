function [mdot_kg_s, h_out_kJkg, regime, Tsat_C, dp_valve_MPa, x_out] = ...
    valve_boiler_to_header(p_up_MPa, p_dn_MPa, opening, x_boiler)
% Valve sizing for steam using Kv with expansion factor (IEC-style approximation).
% Inputs:
%   p_up_MPa   upstream pressure [MPa]
%   p_dn_MPa   downstream pressure [MPa]
%   opening    valve opening fraction [0..1]
%   x_boiler   upstream steam quality [-] (0..1)
% Outputs:
%   mdot_kg_s     mass flow [kg/s]
%   h_out_kJkg    outlet enthalpy [kJ/kg] (isenthalpic throttling)
%   regime        0=no flow, 1=choked, 2=subcritical
%   Tsat_C        upstream saturation temperature [C]
%   dp_valve_MPa  pressure drop across valve [MPa]
%   x_out         downstream quality (recalculated at p_dn)

% ===== Valve & fluid constants =====
Kv_max       = 0.85*1.8e2;   % [m^3/h] initial guess, tune this
Rangeability = 50;      % equal-percentage
Kv_min       = Kv_max / Rangeability;
FL           = 0.90;    % pressure recovery factor
xt           = 0.50;    % choked pressure ratio factor for steam
Fk           = 0.96;    % expansion factor constant for steam

% ===== Outputs init =====
mdot_kg_s    = 0.0;
h_out_kJkg   = 0.0;
regime       = 0.0;
Tsat_C       = 0.0;
dp_valve_MPa = 0.0;
x_out        = x_boiler;

% ===== Input sanitization =====
opening  = max(0.0, min(1.0, double(opening)));
p_up_MPa = double(p_up_MPa);
p_dn_MPa = double(p_dn_MPa);
x_boiler = max(0.0, min(1.0, double(x_boiler)));

% ===== No-flow check =====
if (p_dn_MPa >= p_up_MPa) || (opening <= 0.0)
    p_up_bar   = 10.0 * max(p_up_MPa, 0.0);
    hf = xs_p('hL_p', p_up_bar);
    hg = xs_p('hV_p', p_up_bar);
    h_out_kJkg = (1 - x_boiler)*hf + x_boiler*hg;
    Tsat_C     = xs_p('Tsat_p', p_up_bar);
    regime     = 0.0;
    dp_valve_MPa = 0.0;
    return;
end

% ===== Upstream properties =====
p_up_bar = 10.0 * max(p_up_MPa, 0.0);
Tsat_C   = xs_p('Tsat_p', p_up_bar);
hf       = xs_p('hL_p', p_up_bar);
hg       = xs_p('hV_p', p_up_bar);
vf       = xs_p('vL_p', p_up_bar);
vg       = xs_p('vV_p', p_up_bar);

% Mixture enthalpy and density at boiler outlet
h_in_kJkg = (1 - x_boiler)*hf + x_boiler*hg;
v_mix     = (1 - x_boiler)*vf + x_boiler*vg;
rho_up    = 1.0 / max(v_mix, 1e-9);

% ===== Equal-percentage Kv(opening) =====
Kv_open = min(Kv_max, Kv_min * (Rangeability ^ opening));

% ===== Pressure ratios and absolute drop =====
p_dn_bar = 10.0 * max(p_dn_MPa, 0.0);
dp_bar   = max(p_up_bar - p_dn_bar, 0.0);
x_ratio  = dp_bar / max(p_up_bar, 1e-6);

% ===== Expansion factor and regime =====
x_choke = xt * (FL^2);     % critical pressure ratio
if x_ratio > x_choke
    regime   = 1.0;        % choked
    x_eff    = x_choke;    % use choked ratio for capacity
else
    regime   = 2.0;        % subcritical
    x_eff    = x_ratio;
end

% Effective drop at valve for capacity (bar) and expansion factor
dp_eff_bar = x_eff * p_up_bar;
Y = max(1.0 - x_eff / (3.0 * Fk), 0.0);

% ===== Volumetric flow using Kv at effective Δp (bar) =====
rho_rel = max(rho_up / 1000.0, 1e-6);
Q_m3_h  = Kv_open * Y * sqrt( max(dp_eff_bar, 0.0) / rho_rel );

% ===== Mass flow =====
mdot_kg_s = (Q_m3_h * rho_up) / 3600.0;

% ===== Isenthalpic throttling =====
h_out_kJkg = h_in_kJkg;

% ===== Pressure drop across valve =====
dp_valve_MPa = (p_up_bar - p_dn_bar) / 10.0;  % bar → MPa

% ===== Recalculate downstream quality (optional) =====
p_out_bar = p_dn_bar;
hf_out = xs_p('hL_p', p_out_bar);
hg_out = xs_p('hV_p', p_out_bar);
if h_out_kJkg <= hf_out
    x_out = 0.0; % subcooled
elseif h_out_kJkg >= hg_out
    x_out = 1.0; % dry or superheated
else
    x_out = (h_out_kJkg - hf_out) / (hg_out - hf_out);
end

end

% ===== Local typed wrapper to avoid mxArray propagation =====
function val = xs_p(token, pbar)
coder.extrinsic('XSteam');
val = 0.0;
tmp = feval('XSteam', token, pbar);
if ~isempty(tmp), val = double(tmp); end
if ~isfinite(val), val = 0.0; end
end