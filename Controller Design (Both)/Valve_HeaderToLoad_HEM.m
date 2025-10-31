function [mdot_out_kg_s, h_out_kJkg] = Valve_HeaderToLoad_HEM( ...
    p_up_MPa, x_hdr, P_down_MPa, opening_out)
% Steam_Outflow_Rate
% -------------------------------------------------------------------------
% Homogeneous Equilibrium Model (HEM) outlet (header -> load).
% mdot = Cd * A_eff(u) * sqrt(2 * rho_mix(p,x) * max(ΔP, 0))
% Returns mass flow rate and outlet enthalpy.
% -------------------------------------------------------------------------

% ================= User-tunable valve constants ==========================
CD_OUT        = 0.70;     % [-]  discharge coefficient
AOUT_MAX_M2   = 0.015;    % [m^2] max effective area at u=1 (tune)
R_EQP         = 50;       % [-]  equal-percentage R
USE_EQP_CHAR  = true;     % true: equal-percentage; false: linear A=u*Amax
USE_SOFTPLUS  = false;    % false: hard max; true: stable softplus
DP_SOFT_SCALE = 1e5;      % [Pa] softplus scale if enabled (stable form)
% ========================================================================

coder.extrinsic('XSteam');

% Cast & clamp
p_up_MPa   = double(p_up_MPa);
P_down_MPa = double(P_down_MPa);
x          = max(0.0, min(1.0, double(x_hdr)));
u          = opening_out;

% Effective area
if USE_EQP_CHAR
    f = (R_EQP.^u - 1) / (R_EQP - 1);
else
    f = u;
end
Aeff = max(0.0, AOUT_MAX_M2 * f);

% ΔP in Pa
P_up = max(p_up_MPa,  0) * 1e6;
P_dn = max(P_down_MPa, 0) * 1e6;
if USE_SOFTPLUS
    dP = softplus_stable(P_up - P_dn, DP_SOFT_SCALE);
else
    dP = max(P_up - P_dn, 0.0);
end

if (Aeff <= 0) || (dP <= 0)
    mdot_out_kg_s = 0.0;
    h_out_kJkg = 0.0;
    return;
end

% Mixture properties at header (bar for XSteam)
p_bar = 10.0 * p_up_MPa;
hf = xs('hL_p', p_bar);            % [kJ/kg] saturated liquid enthalpy
hg = xs('hV_p', p_bar);            % [kJ/kg] saturated vapor enthalpy
vf = xs('vL_p', p_bar);            % [m^3/kg]
vg = xs('vV_p', p_bar);            % [m^3/kg]

% Mixture density and enthalpy
v  = max((1.0 - x)*vf + x*vg, 1e-9);     % [m^3/kg]
rho = 1.0 / v;                           % [kg/m^3]
h_out_kJkg = (1.0 - x)*hf + x*hg;        % [kJ/kg]

% HEM orifice flow
mdot_out_kg_s = CD_OUT * Aeff * sqrt( max(2.0 * rho * dP, 0.0) );

end

% -------------------- helpers --------------------
function y = softplus_stable(z, beta)
% Numerically stable softplus: beta*log(1+exp(z/beta))
if z > 0
    y = z + beta * log1p(exp(-z/beta));
else
    y = beta * log1p(exp(z/beta));
end
end

function val = xs(token, pbar)
coder.extrinsic('XSteam');
val = 0.0;
val = double(feval('XSteam', token, pbar));
end
