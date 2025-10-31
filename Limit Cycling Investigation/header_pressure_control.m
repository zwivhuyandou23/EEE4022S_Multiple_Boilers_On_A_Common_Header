% ============================================================
% Header PI tuning via IMC (no load-sharing weights)
% Assumes A,B,C,D exist (linearized model).
% ============================================================

run linearize_system_for_header_pressre_control.m;

sys_full = ss(A,B,C,D);

% Header pressure per each boiler firing rate input
G_Q1 = sys_full(1,1);    % p_hdr per Q1
G_Q2 = sys_full(1,5);    % p_hdr per Q2

% Gains [MPa/MW]
K1 = dcgain(G_Q1);
K2 = dcgain(G_Q2);

% Effective header gain (single demand broadcast to both boilers)
K_eff = K1 + K2;

% Dominant time constant from slowest mode
eigvals = eig(A);
[~,idx_max] = max(real(eigvals));     % closest to zero (negative)
tau = -1/real(eigvals(idx_max));      % [s]

% ============================================================
% Desired closed-loop performance
% ============================================================
Ts_target =  10;    % desired settling time [s]
theta_est = 0.0;    % dead time estimate [s], set from data if known

% Choose lambda (IMC filter parameter)
lambda_req = max(Ts_target/4 - theta_est, 0.1);
lambda_min = 0.25*tau;
lambda     = max(lambda_req, lambda_min);

% ============================================================
% IMC PI gains (FOPDT surrogate)
% ============================================================
Kp = tau / ( K_eff * ( lambda + theta_est ) );
Ti = tau + theta_est;
Ki = Kp / Ti;

ts_pred = 4 * ( lambda + theta_est );

% ============================================================
% Report results
% ============================================================
fprintf('\n=== Header Pressure PI (IMC) ===\n');
fprintf('K1=%.6f, K2=%.6f [MPa/MW], K_eff=%.6f\n', K1, K2, K_eff);
fprintf('tau=%.3f s, theta=%.3f s, lambda=%.3f s\n', tau, theta_est, lambda);
fprintf('Predicted settling ts ≈ %.2f s\n', ts_pred);
fprintf('PI gains: Kp=%.6f, Ki=%.6f  (Ti=%.3f s)\n', Kp, Ki, Ti);

% ============================================================
% Optional: more aggressive alternative
% ============================================================
lambda_fast = max(Ts_target/6 - theta_est, 0.1);
lambda_fast = max(lambda_fast, 0.2*tau);

Kp_fast = tau / ( K_eff * ( lambda_fast + theta_est ) );
Ti_fast = Ti;   % keep integral time tied to plant time constant
Ki_fast = Kp_fast / Ti_fast;
ts_fast = 4 * ( lambda_fast + theta_est );

fprintf('\n-- Aggressive option (validate against limits/noise) --\n');
fprintf('lambda_fast=%.3f s -> ts ≈ %.2f s\n', lambda_fast, ts_fast);
fprintf('PI_fast: Kp=%.6f, Ki=%.6f  (Ti=%.3f s)\n', Kp_fast, Ki_fast, Ti_fast);

% ============================================================
% Optional: simulate closed-loop step response
% ============================================================
s = tf('s');
G_approx = K_eff / (tau*s + 1);
C = pid(Kp, Ki, 0);       % PI controller
CL = feedback(C*G_approx, 1);

figure;
step(CL);
grid on;
title('Closed-loop Step Response (IMC PI)');