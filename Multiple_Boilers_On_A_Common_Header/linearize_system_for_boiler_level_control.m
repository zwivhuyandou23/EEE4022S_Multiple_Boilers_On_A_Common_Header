% === Single Boiler: Improved Tracking & Disturbance Rejection ===
clc; fprintf('=== Single Boiler: Improved PI + Feedforward ===\n');

% Operating point and dynamics (as in your script)
x_star = [57.21; 8.801; 0.08296; 4.900];
u_star = [59.7; 41.0; 41.0; 1290.0];

f_dyn = @(x,u) full_dyn_one_boiler_direct(x,u);
g_out = @(x,u) [ x(1); x(2) ];  % [Level proxy (=Vwt), p_drum]

dx_state = 1e-5 * max(1, abs(x_star));
du_input = zeros(size(u_star));
du_input(2) = 1e-4 * max(1, abs(u_star(2)));  % perturb qf
du_input(3) = 1e-4 * max(1, abs(u_star(3)));  % perturb q4 (disturbance)

[A,B,C,D] = linearize_mimo_central(f_dyn, g_out, x_star, u_star, dx_state, du_input);
sys = ss(A,B,C,D);

% SISO channels
G_L_qf = sys(1,2);     % level per qf (control channel)
G_L_q4 = sys(1,3);     % level per q4 (disturbance channel)

% DC gains
Kqf = dcgain(G_L_qf);
Kq4 = dcgain(G_L_q4);

% Dominant time constant and nominal dead time
poles = pole(G_L_qf);
tau  = -1 / max(real(poles));
theta = 2.0;                   % adjust if sensor/actuator delay known

% Target settling and lambda
Ts_target = 16;                % aim < 20 s for snappier tracking
lambda_req = max(Ts_target/4 - theta, 0.5);   % ensure positive & not too small
lambda_min = 0.25*tau;         % robustness guardrail
lambda = max(lambda_req, lambda_min);

% IMC PI for FOPDT surrogate (qf -> Level)
Kp = tau / ( Kqf * ( lambda + theta ) );
Ti = tau + theta;
Ki = Kp / Ti;

% Two-DoF setpoint weighting
beta = 0.6;                    % proportional acts on weighted setpoint

% Disturbance feedforward gain (evaporation cancellation)
k_ff = - Kq4 / Kqf;

% Bias at operating point
qf_nom  = u_star(2);
q4_nom  = u_star(3);
qf_bias = qf_nom + k_ff * (0); % keep zero if q4 nominal cancels in plant
                               % or set explicitly if a known offset exists

% Report
ts_pred = 4*(lambda + theta);
fprintf('Kqf=%.6f, Kq4=%.6f, tau=%.3f s, theta=%.2f s\n', Kqf, Kq4, tau, theta);
fprintf('lambda=%.2f s -> predicted ts ≈ %.2f s\n', lambda, ts_pred);
fprintf('PI (Two-DoF): Kp=%.4f, Ki=%.5f (Ti=%.2f s), beta=%.2f\n', Kp, Ki, Ti, beta);
fprintf('Feedforward: k_ff=%.4f  (qf_cmd += k_ff * (q4 - q4_nom))\n', k_ff);
fprintf('Bias: qf_bias=%.3f\n', qf_bias);

% End of main script

function xdot = full_dyn_one_boiler_direct(x,u)
% FULL_DYN_ONE_BOILER_DIRECT
% Single-boiler dynamics with feedwater coupling and simple steam generation.

% States: [Vwt; p_drum; x; Vsd]
Vwt = x(1); p_drum = x(2); xq = x(3); Vsd = x(4);

% Inputs: [Q; qf; q4; hf]
Q   = u(1); qf = u(2); q4 = u(3); hf = u(4); %#ok<NASGU>

% Nominals
Vwt_nom = 57.21; p_drum_nom = 8.801; x_nom = 0.08296; Vsd_nom = 4.900;
Q_nom   = 59.7;  qf_nom     = 41.0;   q4_nom = 41.0;

% Steam flow sensitivity to firing (toy coupling)
s_q = 0.5;
q4_eff = q4 + s_q*(Q - Q_nom);

% Drum pressure dynamics
tau_p = 8.0; kQ = 0.015;
p_drum_dot = - (p_drum - p_drum_nom)/tau_p + kQ*(Q - Q_nom);

% Header/quality proxies (single boiler: stable relaxations)
x_dot   = - (xq  - x_nom)/15;
Vsd_dot = - (Vsd - Vsd_nom)/25;

% Feedwater/evaporation coupling into water inventory
k_fw   = 1.0;     % inflow-to-inventory scale (adjust for units)
k_evap = 1.0;     % steam removal scale (proxy via q4_eff)
q_evap     = k_evap * q4_eff;
q_evap_nom = k_evap * q4_nom;

tau_w = 20.0;
Vwt_dot = - (Vwt - Vwt_nom)/tau_w ...
          + k_fw * (qf - qf_nom) ...
          - (q_evap - q_evap_nom);

xdot = [Vwt_dot; p_drum_dot; x_dot; Vsd_dot];
end

function [A,B,C,D] = linearize_mimo_central(f_dyn, g_out, x0, u0, dx_state, du_input)
% LINEARIZE_MIMO_CENTRAL  Central-difference linearization for MIMO systems

n = numel(x0); m = numel(u0); p = numel(g_out(x0,u0));
A = zeros(n,n); B = zeros(n,m); C = zeros(p,n); D = zeros(p,m);

% State Jacobians
for i = 1:n
    dx = zeros(n,1); dx(i) = max(dx_state(i), eps);
    f_plus  = f_dyn(x0 + dx, u0);
    f_minus = f_dyn(x0 - dx, u0);
    A(:,i)  = (f_plus - f_minus) / (2*dx(i));

    y_plus  = g_out(x0 + dx, u0);
    y_minus = g_out(x0 - dx, u0);
    C(:,i)  = (y_plus - y_minus) / (2*dx(i));
end

% Input Jacobians
for j = 1:m
    if du_input(j) == 0, continue; end
    du = zeros(m,1); du(j) = du_input(j);
    f_plus  = f_dyn(x0, u0 + du);
    f_minus = f_dyn(x0, u0 - du);
    B(:,j)  = (f_plus - f_minus) / (2*du(j));

    y_plus  = g_out(x0, u0 + du);
    y_minus = g_out(x0, u0 - du);
    D(:,j)  = (y_plus - y_minus) / (2*du(j));
end
end

function analyze_linearization(A,B,C,D)
% ANALYZE_LINEARIZATION  Print matrices, controllability/observability ranks,
% stability, dominant poles, and DC gain matrix.

fprintf('\n=== Linearization Analysis ===\n');

% Matrices
fprintf('A (%dx%d):\n', size(A,1), size(A,2)); disp(A);
fprintf('B (%dx%d):\n', size(B,1), size(B,2)); disp(B);
fprintf('C (%dx%d):\n', size(C,1), size(C,2)); disp(C);
fprintf('D (%dx%d):\n', size(D,1), size(D,2)); disp(D);

% Controllability / Observability
Co = ctrb(A,B); Ob = obsv(A,C);
fprintf('Controllability rank: %d of %d\n', rank(Co), size(A,1));
fprintf('Observability   rank: %d of %d\n', rank(Ob), size(A,1));

% Stability and dominant pole
eigvals = eig(A);
[~,idx_max] = max(real(eigvals)); p_dom = eigvals(idx_max);
isStable = all(real(eigvals) < 0);
fprintf('Stability: %s (max Re(λ)=%.6f)\n', ternary(isStable,'Stable','Unstable'), max(real(eigvals)));
if real(p_dom) < 0
    tau_dom = -1/real(p_dom);
    fprintf('Dominant pole λ = %.6f, tau ≈ %.3f s\n', p_dom, tau_dom);
else
    fprintf('Dominant pole λ = %.6f (unstable)\n', p_dom);
end
fprintf('Eigenvalues:\n'); disp(eigvals);

% DC gain matrix
sys_full = ss(A,B,C,D);
K_dc = dcgain(sys_full);
fprintf('DC gain matrix (outputs × inputs):\n'); disp(K_dc);
end

function s = ternary(cond, a, b)
% TERNARY  Inline conditional
if cond, s = a; else, s = b; end
end