function [dpdt_MPa, dxdt] = Common_Header_SegmentBC( ...
    p_hdr_MPa, x, m_in_vec, h_in_vec, m_out_kg_s)
% COMMON_HEADER_SEGMENTBC (Simulink-safe, improved)
% -------------------------------------------------------------------------
% Wet-steam common header dynamics at saturation (constant volume).
% States (integrated externally):
%   p_hdr_MPa [MPa] : header average pressure
%   x         [-]   : steam quality (mass fraction vapor), clamped 0..1
%
% Inputs:
%   m_in_vec [kg/s]  : vector of inlet mass flows (>=0)
%   h_in_vec [kJ/kg] : vector of inlet enthalpies (same length as m_in_vec)
%   m_out_kg_s [kg/s]: outlet mass flow (>=0), enthalpy assumed h_mix (header)
%
% Outputs:
%   dpdt_MPa [MPa/s] : pressure time derivative
%   dxdt     [1/s]   : quality time derivative
%
% Assumptions:
% - Saturation mixture in the header, properties from XSteam via typed wrapper.
% - Energy balance includes -V * dp work (1 bar·m^3/s = 100 kJ/s).
% - Outlet enthalpy equals header mixture enthalpy (isenthalpic draw from header).
% -------------------------------------------------------------------------

% ================= User-tunable geometry =================================
D_hdr   = 1.20;       % [m] header diameter (inner)
L_hdr   = 50.0;       % [m] header length
% Roughness is not used here since friction DP is not returned; keep for future
eps_hdr = 4.5e-5;     % [m] roughness (steel) (placeholder)

% Derived geometry
A_hdr       = pi * D_hdr^2 / 4;   % [m^2] cross-sectional area
V_HEADER_M3 = A_hdr * L_hdr;      % [m^3] control volume for header

% ================= Numerical constants ===================================
DP_STEP_REL   = 1e-3;
DP_STEP_MIN   = 1e-4;
REG_EPS_BASE  = 1e-10;
P_MIN_BAR     = 1.0;
P_MAX_BAR     = 200.0;
MU_NOM        = 1e-5;  %#ok<NASGU> % Pa·s, kept for potential future friction modeling
% ========================================================================

coder.extrinsic('XSteam');

% ---- Types & clamps
p_hdr_MPa   = double(p_hdr_MPa);
x           = max(0.0, min(1.0, double(x)));
m_out_kg_s  = max(0.0, double(m_out_kg_s));

% ---- Fixed-size handling for Simulink Coder (no dynamic allocation)
nin = min(numel(m_in_vec), numel(h_in_vec));
if nin <= 0
    % No inlets -> pure blowdown; keep types
    dpdt_MPa = 0.0;
    dxdt     = 0.0;
    return;
end
m_in = zeros(nin,1);
h_in = zeros(nin,1);
for i = 1:nin
    m_in(i) = max(0.0, double(m_in_vec(i)));
    h_in(i) = double(h_in_vec(i));
end

% ---- Early out if volume invalid
if ~(V_HEADER_M3 > 0)
    dpdt_MPa = 0.0; dxdt = 0.0; return;
end

% ---- Properties pressure in bar
p_bar_raw = 10.0 * max(p_hdr_MPa, 0.0);
p_bar     = min(max(p_bar_raw, P_MIN_BAR), P_MAX_BAR);

% ---- Saturation properties (typed wrapper)
h_f = xs('hL_p', p_bar);  % kJ/kg
h_g = xs('hV_p', p_bar);  % kJ/kg
v_f = xs('vL_p', p_bar);  % m^3/kg
v_g = xs('vV_p', p_bar);  % m^3/kg

% ---- Mixture properties
h_fg  = max(h_g - h_f, 1e-8);
h_mix = h_f + x * h_fg;
v_mix = max((1.0 - x) * v_f + x * v_g, 1e-9);
rho   = 1.0 / v_mix;                  % kg/m^3
Mhdr  = rho * V_HEADER_M3;            % kg

% ---- Finite-difference derivatives wrt p (bar)
dpb  = max(DP_STEP_MIN, DP_STEP_REL * max(p_bar, 1.0));
p_lo = max(p_bar - dpb, P_MIN_BAR);
p_hi = min(p_bar + dpb, P_MAX_BAR);
if p_hi <= p_lo
    if p_bar <= P_MIN_BAR + 10*DP_STEP_MIN
        p_hi = min(p_bar + DP_STEP_MIN, P_MAX_BAR);
        p_lo = p_bar;
    else
        p_lo = max(p_bar - DP_STEP_MIN, P_MIN_BAR);
        p_hi = p_bar;
    end
end
denom = max(p_hi - p_lo, 1e-8);

% Saturation properties at perturbed pressures
h_f_m = xs('hL_p', p_lo);  h_f_p = xs('hL_p', p_hi);
h_g_m = xs('hV_p', p_lo);  h_g_p = xs('hV_p', p_hi);
v_f_m = xs('vL_p', p_lo);  v_f_p = xs('vL_p', p_hi);
v_g_m = xs('vV_p', p_lo);  v_g_p = xs('vV_p', p_hi);

% Finite-difference derivatives w.r.t pressure
dhL_dpbar = (h_f_p - h_f_m) / denom;
dhV_dpbar = (h_g_p - h_g_m) / denom;
dvL_dpbar = (v_f_p - v_f_m) / denom;
dvV_dpbar = (v_g_p - v_g_m) / denom;

% Partials at mixture
dhdp_x_bar   = dhL_dpbar + x * (dhV_dpbar - dhL_dpbar);
dhdx_p       = h_fg;
dvdp_mix     = (1.0 - x) * dvL_dpbar + x * dvV_dpbar;
drhodp_x_bar = - dvdp_mix / (v_mix^2);
drhodx_p     = - (v_g - v_f) / (v_mix^2);

% ---- Mass & energy rates
Min  = sum(m_in);
Mout = m_out_kg_s;
mdot = Min - Mout;                 % net mass accumulation [kg/s]
Ein  = sum(m_in .* h_in);          % kJ/s (since h is kJ/kg)
Eout = Mout * h_mix;               % kJ/s (header mixture enthalpy assumption)

% ---- Linear system for [dpdt_bar; dxdt]
% Balance:
%   dM/dt = mdot
%   d(M*h_mix)/dt = Ein - Eout - V * (dp/dt_bar) * 100
% Expand energy derivative:
%   h_mix * dM/dt + Mhdr * (dhdp_x_bar * dpdt_bar + dhdx_p * dxdt) = Ein - Eout - V*100*dpdt_bar
%
% Rearranged to A * [dpdt_bar; dxdt] = B
A11 = V_HEADER_M3 * drhodp_x_bar;                         % from dM/dt = V * drho/dp * dpdt_bar + V*drho/dx * dxdt, but we keep mass eq in B1
A12 = V_HEADER_M3 * drhodx_p;

% Energy equation coefficients
A21 = Mhdr * dhdp_x_bar - V_HEADER_M3 * 100.0;           % 100 kJ/s per bar·m^3/s
A22 = Mhdr * dhdx_p;

% Right-hand sides
B1  = mdot;                                              % mass balance
B2  = Ein - Eout - h_mix * mdot;                         % energy balance after moving h_mix*dM/dt

% Solve 2x2 (regularized)
detA = A11*A22 - A12*A21;
REG  = REG_EPS_BASE;
if ~isfinite(detA) || abs(detA) < 1e-9
    REG = 1e-6;
end
A = [A11, A12; A21, A22];
B = [B1;  B2];

sol = (A + REG*eye(2)) \ B;

dpdt_bar = sol(1);
dxdt     = sol(2);

% ---- Project dxdt at the bounds (no leave the [0,1] interval)
if x <= 0.0 + 1e-6 && dxdt < 0.0
    dxdt = 0.0;
elseif x >= 1.0 - 1e-6 && dxdt > 0.0
    dxdt = 0.0;
end

% Final safety
if ~isfinite(dpdt_bar), dpdt_bar = 0.0; end
if ~isfinite(dxdt),     dxdt     = 0.0; end

% Convert to MPa/s from bar/s
dpdt_MPa = dpdt_bar / 10.0;
end

% ================= Local typed helper for XSteam =========================
function val = xs(token, pbar)
coder.extrinsic('XSteam');
val = 0.0;
tmp = feval('XSteam', token, pbar);
if ~isempty(tmp)
    val = double(tmp);
end
if ~isfinite(val)
    val = 0.0;
end
end