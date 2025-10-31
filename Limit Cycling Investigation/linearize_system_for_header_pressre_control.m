%% linearize_system_all.m
% System-wide linearization of two boilers feeding a common header (MIMO).

clc;
%% === Steady operating point (from user) ===
% States: [Vwt1; p_drum1; x1; Vsd1; Vwt2; p_drum2; x2; Vsd2; p_hdr; x_hdr]
x_star = [57.24; 8.807; 0.08298; 4.901; ...
          57.30; 8.693; 0.06649; 4.901; ...
          8.513; 0.08218];

% Inputs: [Q1; qf1; q4_1; hf1; Q2; qf2; q4_2; hf2; m_out]
u_star = [59.7; 41.0; 41.0; 1290.0; ...
          49.9; 34.2; 34.2; 1290.0; ...
          71.97];

%% === Plant and output wrappers (direct-input model) ===
f_dyn = @(x,u) full_dyn_two_boiler_header_direct(x,u);
g_out = @(x,u) [x(9); x(2); x(6)];  % [p_hdr; p_drum1; p_drum2]

%% === Perturbation magnitudes (state and input) ===
dx_state = 1e-5 * max(1, abs(x_star));   % per-state perturbations

% Only Q1, qf1, Q2, qf2 get nonzero perturbations; others are fixed (0)
du_input = zeros(size(u_star));
du_input(1) = 1e-4 * max(1, abs(u_star(1)));  % Q1
%du_input(2) = 1e-4 * max(1, abs(u_star(2)));  % qf1
du_input(5) = 1e-4 * max(1, abs(u_star(5)));  % Q2
%du_input(6) = 1e-4 * max(1, abs(u_star(6)));  % qf2
% du_input(3)=0; du_input(4)=0; du_input(7)=0; du_input(8)=0; du_input(9)=0

%% === MIMO central-difference linearization ===
[A,B,C,D] = linearize_mimo_central(f_dyn, g_out, x_star, u_star, dx_state, du_input);

%% === Diagnostics ===
print_lin_diagnostics(A,B,C,D);

Co = ctrb(A,B); Ob = obsv(A,C);
fprintf('Controllability rank: %d of %d\n', rank(Co), size(A,1));
fprintf('Observability   rank: %d of %d\n', rank(Ob), size(A,1));

poles = eig(A); isStable = all(real(poles)<0);
fprintf('Stability: %s (max Re(λ)=%.6f)\n', ternary(isStable,'Stable','Unstable'), max(real(poles)));

[~,idx_max] = max(real(poles)); p_dom = poles(idx_max);[~,idx_max] = max(real(poles));
p_dom = poles(idx_max);
tau_dom = -1 / real(p_dom);  % valid when real(p_dom) < 0
fprintf('Dominant pole λ=%.6f, tau≈%.3f s\n', p_dom, tau_dom);

sys = ss(A,B,C,D);
K_dc = dcgain(sys);   % full MIMO DC gain matrix (outputs × inputs)
disp('DC gain matrix (outputs × inputs):'); disp(K_dc);

%% === Helper: inline ternary for printing ===
function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function xdot = full_dyn_two_boiler_header_direct(x,u)
% FULL_DYN_TWO_BOILER_HEADER_DIRECT
% Minimal, stable dynamics for linearization with direct physical inputs.
% State vector (10x1): [Vwt1; p_drum1; x1; Vsd1; Vwt2; p_drum2; x2; Vsd2; p_hdr; x_hdr]
% Input vector (9x1):  [Q1; qf1; q4_1; hf1; Q2; qf2; q4_2; hf2; m_out]

% --- Unpack states
Vwt1=x(1); p_drum1=x(2); x1=x(3); Vsd1=x(4);
Vwt2=x(5); p_drum2=x(6); x2=x(7); Vsd2=x(8);
p_hdr=x(9); x_hdr=x(10);

% --- Unpack inputs
Q1=u(1); qf1=u(2); q4_1=u(3); hf1=u(4);
Q2=u(5); qf2=u(6); q4_2=u(7); hf2=u(8);
m_out=u(9);

% --- Nominal values
Vwt1_nom=57.21; p_drum1_nom=8.801; x1_nom=0.08296; Vsd1_nom=4.900;
Vwt2_nom=57.43; p_drum2_nom=8.700; x2_nom=0.06648; Vsd2_nom=4.900;
p_hdr_nom=8.521; x_hdr_nom=0.08235;
Q1_nom=59.7; Q2_nom=49.9; q4_1_nom=41.0; q4_2_nom=34.2;

% --- Steam flow sensitivities to firing (toy first-order coupling)
s_q1=0.5; s_q2=0.6;
q4_1_eff = q4_1 + s_q1*(Q1 - Q1_nom);
q4_2_eff = q4_2 + s_q2*(Q2 - Q2_nom);

% --- Header pressure dynamics (mass balance + compression)
p_hdr_ref = 0.5*(p_drum1 + p_drum2);
k_p_mass  = 1e-3; tau_hdr = 5.0;
p_hdr_dot = k_p_mass*((q4_1_eff + q4_2_eff) - m_out) - (p_hdr - p_hdr_ref)/tau_hdr;

% --- Header quality proxy
x_hdr_dot = -(x_hdr - x_hdr_nom)/10.0;

% --- Drum pressure dynamics: first-order response to firing and relaxation
tau_p1=8.0; tau_p2=7.0; kQ1=0.015; kQ2=0.018;
p_drum1_dot = - (p_drum1 - p_drum1_nom)/tau_p1 + kQ1*(Q1 - Q1_nom);
p_drum2_dot = - (p_drum2 - p_drum2_nom)/tau_p2 + kQ2*(Q2 - Q2_nom);

% --- Boiler inventory/quality/steam-dome volumes: stable relaxation to nominal
Vwt1_dot   = - (Vwt1 - Vwt1_nom)/20;
x1_dot     = - (x1   - x1_nom)/15;
Vsd1_dot   = - (Vsd1 - Vsd1_nom)/25;

Vwt2_dot   = - (Vwt2 - Vwt2_nom)/20;
x2_dot     = - (x2   - x2_nom)/15;
Vsd2_dot   = - (Vsd2 - Vsd2_nom)/25;

xdot = [Vwt1_dot; p_drum1_dot; x1_dot; Vsd1_dot; ...
        Vwt2_dot; p_drum2_dot; x2_dot; Vsd2_dot; ...
        p_hdr_dot; x_hdr_dot];
end

function [A,B,C,D] = linearize_mimo_central(f_dyn, g_out, x0, u0, dx_state, du_input)
% LINEARIZE_MIMO_CENTRAL  Central-difference linearization for MIMO systems
f0 = f_dyn(x0,u0); %#ok<NASGU>
y0 = g_out(x0,u0); %#ok<NASGU>

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

% Input Jacobians (respect per-input perturbations)
for j = 1:m
    if du_input(j) == 0
        B(:,j) = 0;
        D(:,j) = 0;
        continue;
    end
    du = zeros(m,1); du(j) = du_input(j);
    f_plus  = f_dyn(x0, u0 + du);
    f_minus = f_dyn(x0, u0 - du);
    B(:,j)  = (f_plus - f_minus) / (2*du(j));

    y_plus  = g_out(x0, u0 + du);
    y_minus = g_out(x0, u0 - du);
    D(:,j)  = (y_plus - y_minus) / (2*du(j));
end
end

function print_lin_diagnostics(A,B,C,D)
fprintf('\n=== Linearization Diagnostics ===\n');
fprintf('A: %dx%d, B: %dx%d\n', size(A,1), size(A,2), size(B,1), size(B,2));
fprintf('C: %dx%d, D: %dx%d\n', size(C,1), size(C,2), size(D,1), size(D,2));

lambda = eig(A);
[~, order] = sort(real(lambda), 'ascend');
lambda = lambda(order);

fprintf('Eigenvalues (real, imag):\n');
for k = 1:numel(lambda)
    fprintf('  λ%-2d = %+ .6f  %+ .6f i\n', k, real(lambda(k)), imag(lambda(k)));
end

rp = real(lambda).';
ip = imag(lambda).';
fprintf('Real parts:\n'); disp(rp);
fprintf('Imag parts:\n'); disp(ip);

fprintf('\nA matrix:\n'); disp(A);
fprintf('B matrix:\n'); disp(B);
fprintf('C matrix:\n'); disp(C);
fprintf('D matrix:\n'); disp(D);
end