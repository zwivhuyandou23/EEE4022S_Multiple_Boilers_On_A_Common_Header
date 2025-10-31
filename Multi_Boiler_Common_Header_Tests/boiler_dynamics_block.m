function [dxdt, additional_outputs] = boiler_dynamics_block(x, u, b, k, Vr, Adc, m_metal, m_riser, tau_s)
% Åström & Bell (2000) 4-state drum-boiler model (extended diagnostics)
% States:  x = [Vwt; p_MPa; a3; Vsd]
% Inputs:  u = [Q_MW; qf; q4; hf]
% Outputs: dxdt = [dVwt; dp; da3; dVsd]
%          additional_outputs = [qdc; q3; CR; res_eff; e33; e32dp_raw]

coder.extrinsic('get_boiler_calcs');

% --- Unpack state and input ---
Vwt   = x(1); 
p_MPa = x(2); 
a3    = x(3); 
Vsd   = x(4);

Q_MW  = u(1); 
qf    = u(2); 
q4    = u(3); 
hf    = u(4);

% --- Constants ---
Vd   = 40; 
Vdc  = 11; 
Vsd0 = 4.9;
Cp_m = 0.5; 
g    = 9.81; 
EPS  = 1e-12;
kp = 1;

dxdt               = zeros(4,1);
additional_outputs = zeros(6,1);

% --- Get steam properties and void fraction ---
props = struct('rhoL',0,'rhoV',0,'hL',0,'hV',0,'vL',0,'vV',0,'Tsat',0, ...
               'drhoL_dpbar',0,'drhoV_dpbar',0,'dhL_dpbar',0,'dhV_dpbar',0, ...
               'dvL_dpbar',0,'dvV_dpbar',0,'dTs_dpbar',0);
av = 0; dav_dp = 0; dav_da3 = 0;
[props, av, dav_dp, dav_da3] = get_boiler_calcs(p_MPa, a3);

% --- Convert derivatives from per bar to per MPa ---
rho_w   = props.rhoL;  drho_w_dp = props.drhoL_dpbar * 10;
rho_s   = props.rhoV;  drho_s_dp = props.drhoV_dpbar * 10;
h_w     = props.hL;    dh_w_dp  = props.dhL_dpbar    * 10;
h_s     = props.hV;    dh_s_dp  = props.dhV_dpbar    * 10;
dTs_dp  = props.dTs_dpbar       * 10;

% --- Guards ---
a3      = min(max(a3, 1e-6), 0.30);
av      = min(max(av, 0.0), 1.0);
dav_dp  = min(max(dav_dp, -1e6), 1e6);
dav_da3 = min(max(dav_da3, 1e-8), 1e6);

% --- Thermo aggregates ---
h_fg = max(h_s - h_w, 1e-6);
Vt   = Vd + Vr + Vdc;
Vs   = Vt - Vwt;
Q    = Q_MW * 1000;  % MW → kJ/s

% --- Circulation (Eq. 15) ---
circ_num = 2 * rho_w * Adc^2 * max(rho_w - rho_s, 0) * g * av * Vr;
qdc      = sqrt(max(circ_num / max(k, EPS), 0));

% ================== Coefficients (Eq. 25) ==================
% --- 2x2 block for [dVwt; dp] ---
e11 = rho_w - rho_s;
e12 = Vwt * drho_w_dp + Vs * drho_s_dp;

e21 = rho_w * h_w - rho_s * h_s;
e22 = Vwt * (h_w * drho_w_dp + rho_w * dh_w_dp) + ...
      Vs  * (h_s * drho_s_dp + rho_s * dh_s_dp) - ...
      Vt  * kp + m_metal * Cp_m * dTs_dp;

% --- Riser quality dynamics ---
e33 = h_fg * Vr * ((1 - a3) * rho_s + a3 * rho_w) * dav_da3;

e32 = (rho_w * dh_w_dp - a3 * h_fg * drho_w_dp) * (1 - av) * Vr + ...
      (rho_s * dh_s_dp + (1 - a3) * h_fg * drho_s_dp) * av * Vr + ...
      (rho_s - rho_w) * h_fg * Vr * dav_dp + ...
      m_riser * Cp_m * dTs_dp;

% --- Drum steam volume dynamics ---
e44 = rho_s;

e42 = Vsd * drho_s_dp + ...
      (1 / h_fg) * (rho_s * Vsd * dh_s_dp + rho_w * Vwt * dh_w_dp - Vsd) - ...
      Vwt + m_riser * Cp_m * dTs_dp + ...
      a3 * (1 + b) * Vr * (av * drho_s_dp + (1 - av) * drho_w_dp + (rho_s - rho_w) * dav_dp);

e43 = a3 * (1 + b) * (rho_s - rho_w) * Vr * dav_da3;

% ================== Solve sequentially ==================
% --- (i) [dVwt; dp] ---
A12   = [e11, e12; e21, e22];
rhs12 = [qf - q4; Q + qf * hf - q4 * h_s];
if rcond(A12) < 1e-14
    dVwt = 0; dp = 0;
else
    sol12 = A12 \ rhs12;
    dVwt  = sol12(1);
    dp    = sol12(2);
end

% --- (ii) da3 ---
rhs3 = Q - a3 * h_fg * qdc;
if abs(e33) < 1e-12
    da3 = 0;
else
    da3 = (rhs3 - e32 * dp) / e33;
end

% --- (iii) dVsd ---
rhs4 = (rho_s / max(tau_s, EPS)) * (Vsd0 - Vsd);
if abs(e44) < 1e-12
    dVsd = 0;
else
    dVsd = (rhs4 - e42 * dp - e43 * da3) / e44;
end

dxdt = [dVwt; dp; da3; dVsd];

% ================== Diagnostics ==================
% Condensation flow (Eq. 24)
q_ct = (1 / max(h_fg, EPS)) * ( ...
    Vwt * (h_w * drho_w_dp + rho_w * dh_w_dp) + ...
    Vs  * (h_s * drho_s_dp + rho_s * dh_s_dp) + ...
    m_metal * Cp_m * dTs_dp ) * dp;

% Riser steam flow to drum (Eq. 22)
q3 = qdc ...
   - Vr * (av * drho_s_dp + (1 - av) * drho_w_dp + (rho_s - rho_w) * dav_dp) * dp ...
   - (rho_w - rho_s) * Vr * dav_da3 * da3;

% Circulation ratio
CR = qdc / max(q4, EPS);

% "Residence efficiency" example metric
res_eff = (q3 - q_ct) / max(q3, EPS);

% Raw e32*dp term for debugging
e32dp_raw = e32 * dp;

% --- Extended outputs ---
additional_outputs = [qdc; q3; CR; res_eff; e33; e32dp_raw];
end