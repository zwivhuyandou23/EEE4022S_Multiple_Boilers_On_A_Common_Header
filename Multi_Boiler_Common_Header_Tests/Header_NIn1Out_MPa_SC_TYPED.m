function [dpdt_MPa, dxdt, p_out_MPa, dp_fric_MPa, h_mix_out] = ...
  Header_NIn1Out_MPa_SC_TYPED(p_hdr_MPa, x, m_out_kg_s, m_in_vec, h_in_vec, geom, loss)
% LargeHeaderFcn
% ----------------
% Dynamic wet-steam header (saturated, constant-volume) with:
%   - Two-state storage dynamics: pressure p and quality x
%   - Mass and energy balances with XSteam properties
%   - Darcy–Weisbach frictional drop (Reynolds-dependent) + minor losses
%   - Effective outlet pressure p_out_MPa for the consumer (load)
%
% This function is designed for a "large trunk header" downstream of a mixing
% node (e.g., Node B). It takes one outflow (to the load) and N inflows
% (combined feeds from upstream branches/boilers).
%
% Inputs
%   p_hdr_MPa   [MPa]    Header pressure (state variable)
%   x           [-]      Header steam quality (state variable, 0..1)
%   m_out_kg_s  [kg/s]   Outflow to the consumer (load) from this header
%   m_in_vec    [kg/s]   Column vector of inflow mass rates to the header
%   h_in_vec    [kJ/kg]  Column vector of inflow enthalpies (aligned to m_in_vec)
%   geom                 Struct for geometry (optional; defaults provided)
%       .D     [m]       Inside diameter of the header (default: 0.40)
%       .L     [m]       Header length (default: 40.0)
%   loss                 Struct for hydraulic loss parameters (optional)
%       .rough_m       [m]    Absolute roughness (steel ~0.045 mm; default: 4.5e-5)
%       .K_minor       [-]    Lumped minor loss coefficient (default: 0.5)
%       .mu_Pa_s       [Pa·s] Dynamic viscosity (constant approximation; default: 1.5e-5)
%
% Outputs
%   dpdt_MPa     [MPa/s]  Pressure rate (storage dynamics)
%   dxdt         [1/s]    Quality rate
%   k_MPa_per_kg [MPa/kg] Disturbance gain (P(s)/D(s) = k/s) w.r.t. m_out
%   p_out_MPa    [MPa]    Effective outlet pressure after friction/minor losses
%   dp_fric_MPa  [MPa]    Total Darcy–Weisbach + minor loss across header to load
%   h_mix_out    [kJ/kg]  Mixture enthalpy at header outlet (used by downstream mixing)
%
% Notes
%   - The header is modeled as a well-mixed saturated control volume.
%   - Thermodynamic properties and derivatives use XSteam at the current pressure.
%   - The hydraulic drop is an algebraic effect over the header length to the load.
%   - Use p_out_MPa as the downstream pressure for connected valves/consumer.
%   - This function is compatible with code generation by guarding XSteam calls.
%
% Example usage
%   geom = struct('D',0.10,'L',40.0);
%   loss = struct('rough_m',4.5e-5,'K_minor',0.5,'mu_Pa_s',1.5e-5);
%   [dpdt, dxdt, k, p_out, dp_fric, h_out] = LargeHeaderFcn(8.0, 0.1, 30.0, [15;20], [3200;3100], geom, loss);

arguments
  p_hdr_MPa   (1,1) double
  x           (1,1) double
  m_out_kg_s  (1,1) double
  m_in_vec    (:,1) double
  h_in_vec    (:,1) double
  geom        struct = struct()
  loss        struct = struct()
end

% ---- Defaults (self-contained, overridable via geom/loss) ----
D_hdr = getfield_or(geom, 'D', 0.25);      % [m]
L_hdr = getfield_or(geom, 'L', 50.0);      % [m]
A_hdr = pi * D_hdr^2 / 4;
V_hdr = A_hdr * L_hdr;

rough_m = getfield_or(loss, 'rough_m', 4.5e-5);   % [m] steel
K_minor = getfield_or(loss, 'K_minor', 0.5);      % [-] tees/elbows lumped
mu_Pa_s = getfield_or(loss, 'mu_Pa_s', 1.5e-5);   % [Pa·s] steam (approx)

% ---- Pre-allocated outputs (for codegen shape stability) ----
dpdt_MPa     = 0.0;
dxdt         = 0.0;
k_MPa_per_kg = 0.0;
p_out_MPa    = max(p_hdr_MPa, 0.0);
dp_fric_MPa  = 0.0;
h_mix_out    = 0.0;

% ---- Input sanitation ----
x          = max(0.0, min(1.0, x));
m_out_kg_s = max(0.0, m_out_kg_s);

nin = min(numel(m_in_vec), numel(h_in_vec));
if nin <= 0 || V_hdr <= 0
  % No inflows or invalid geometry -> remain at defaults
  return;
end

% ---- Copy and clamp inflows ----
m_in = zeros(nin,1);
h_in = zeros(nin,1);
for i = 1:nin
  m_in(i) = max(0.0, m_in_vec(i));
  h_in(i) =        h_in_vec(i);
end

% ---- Thermodynamic properties at current pressure (XSteam) ----
% We treat the header mixture as saturated at p_hdr.
p_bar_raw = 10.0 * max(p_hdr_MPa, 0.0);     % [bar]
p_bar     = min(max(p_bar_raw, 1.0), 200.0);% clamp to XSteam table range

h_f = xs('hL_p', p_bar);    % [kJ/kg] saturated liquid enthalpy
h_g = xs('hV_p', p_bar);    % [kJ/kg] saturated vapor enthalpy
v_f = xs('vL_p', p_bar);    % [m^3/kg] saturated liquid specific volume
v_g = xs('vV_p', p_bar);    % [m^3/kg] saturated vapor specific volume

% ---- Mixture properties ----
h_fg  = max(h_g - h_f, 1e-8);
h_mix = h_f + x * h_fg;
v_mix = max((1 - x)*v_f + x*v_g, 1e-9);
rho   = 1.0 / v_mix;     % [kg/m^3]
M_hdr = rho * V_hdr;     % [kg] mass in header

% ---- Finite-difference property derivatives wrt pressure ----
% For dp/dt coupling, we need ∂ρ/∂p and ∂h/∂p at the current x.
dpb   = max(1e-4, 1e-3 * max(p_bar,1.0));  % [bar] step
p_lo  = max(p_bar - dpb, 1.0);
p_hi  = min(p_bar + dpb, 200.0);
denom = max(p_hi - p_lo, 1e-8);

hL_lo = xs('hL_p', p_lo); hL_hi = xs('hL_p', p_hi);
hV_lo = xs('hV_p', p_lo); hV_hi = xs('hV_p', p_hi);
vL_lo = xs('vL_p', p_lo); vL_hi = xs('vL_p', p_hi);
vV_lo = xs('vV_p', p_lo); vV_hi = xs('vV_p', p_hi);

dhL = (hL_hi - hL_lo) / denom;
dhV = (hV_hi - hV_lo) / denom;
dvL = (vL_hi - vL_lo) / denom;
dvV = (vV_hi - vV_lo) / denom;

dhdp_x_bar   = dhL + x*(dhV - dhL);               % ∂h/∂p at fixed x
dhdx_p       = h_fg;                               % ∂h/∂x at fixed p
dvdp_mix     = (1 - x)*dvL + x*dvV;                % ∂v/∂p (mix)
drhodp_x_bar = - dvdp_mix / (v_mix^2);             % ∂ρ/∂p via ρ = 1/v
drhodx_p     = - (v_g - v_f) / (v_mix^2);          % ∂ρ/∂x at fixed p

% ---- Mass & energy balances ----
Min  = sum(m_in);           % [kg/s]
Mout = m_out_kg_s;          % [kg/s]
mdot = Min - Mout;          % net mass rate [kg/s]

Ein  = sum(m_in .* h_in);   % [kJ/s]
Eout = Mout * h_mix;        % [kJ/s]

% ---- Linear system for storage dynamics ----
% A * [dpdt_bar; dxdt] = B, where p is in [bar] internally for XSteam compatibility.
A11 = V_hdr * drhodp_x_bar;
A12 = V_hdr * drhodx_p;
A21 = M_hdr * dhdp_x_bar - V_hdr * 100.0;  % unit-consistency term (matches your base)
A22 = M_hdr * dhdx_p;

A = [A11, A12; A21, A22];
B = [ mdot;
      Ein - Eout - h_mix*mdot ];

% Regularization for numerical robustness
REG  = 1.0e-10;
detA = A11*A22 - A12*A21;
if ~isfinite(detA) || abs(detA) < 1e-9
  REG = 1.0e-6;
end

sol      = (A + REG*eye(2)) \ B;
dpdt_bar = sol(1);     % [bar/s]
dxdt     = sol(2);     % [1/s]

% Clamp dxdt near bounds
if x <= 1e-6  && dxdt <  0, dxdt = 0; end
if x >= 1-1e-6 && dxdt >  0, dxdt = 0; end

% Convert to MPa/s
dpdt_MPa = dpdt_bar / 10.0;

% ---- Disturbance gain k via central difference in Mout ----
delta = 1e-3;  % [kg/s] perturbation
B2_const = Ein - h_mix * Min;

% +delta
Mout_p     = Mout + delta;
B1_plus    = Min - Mout_p;
sol_plus   = (A + REG*eye(2)) \ [B1_plus; B2_const];
dpdt_plus  = sol_plus(1);   % [bar/s]

% -delta
Mout_m     = Mout - delta;
B1_minus   = Min - Mout_m;
sol_minus  = (A + REG*eye(2)) \ [B1_minus; B2_const];
dpdt_minus = sol_minus(1);  % [bar/s]

sens_bar_per_kg_s = (dpdt_plus - dpdt_minus) / (2*delta);
k_MPa_per_kg      = sens_bar_per_kg_s / 10.0;  % [MPa/kg]

% ---- Darcy–Weisbach friction + minor losses to the load ----
% Use outflow Mout, header area A_hdr, and mixture density rho.
v  = Mout / (rho * A_hdr);                                  % [m/s]
Re = max(rho * v * D_hdr / max(mu_Pa_s,1e-9), 1.0);         % Reynolds #

% Friction factor: laminar vs Swamee–Jain turbulent
if Re < 2300.0
  f_D = 64.0 / max(Re, 1.0e-6);
else
  rr  = max(rough_m / max(D_hdr,1.0e-9), 0.0);              % ε/D
  f_D = ( -2.0 * log10( rr/3.7 + 5.74 / (Re^0.9) ) )^(-2);
  f_D = min(max(f_D, 5e-4), 0.2);                           % sanity clamp
end

dp_fric_Pa  = f_D     * (L_hdr / D_hdr) * (rho * v^2 / 2.0); % [Pa]
dp_minor_Pa = K_minor *               1.0 * (rho * v^2 / 2.0);% [Pa]
dp_fric_MPa = (dp_fric_Pa + dp_minor_Pa) / 1.0e6;            % [MPa]

% Effective outlet pressure at the consumer (load)
p_out_MPa = max(p_hdr_MPa - dp_fric_MPa, 0.0);

% Header outlet enthalpy (well-mixed saturated at current state)
h_mix_out = h_mix;

end

% =========================
% Helper utilities (local)
% =========================
function val = getfield_or(s, f, d)
  if ~isempty(s) && isstruct(s) && isfield(s, f)
    val = s.(f);
  else
    val = d;
  end
end

function val = xs(token, pbar)
% Safe wrapper for XSteam(token, pbar) with codegen guard.
  coder.extrinsic('XSteam');
  val = 0.0;
  tmp = feval('XSteam', token, pbar);
  if ~isempty(tmp), val = double(tmp); end
  if ~isfinite(val), val = 0.0; end
end