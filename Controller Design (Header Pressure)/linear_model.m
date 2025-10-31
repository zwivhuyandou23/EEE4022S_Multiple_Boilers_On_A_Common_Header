% Make sure boiler_dynamics_block.m and linearize_boiler.m are in your path
% Initialise both boilers
run Init_two_Boilers.m   % should define x0_1,u0_1 and x0_2,u0_2

function [A,B,E,C,eqns,diag] = linearise_two_boiler_header_full(useSeriesValvePipe)
% Linearise the two-boiler + dynamic header plant and print explicit equations.
% Inputs:
%   useSeriesValvePipe : true  -> use SeriesValvePipe_HEM_DW_SC for valve inflow sensitivities
%                        false -> use valve_boiler_to_header for valve inflow sensitivities
%
% Outputs:
%   A,B,E,C : global continuous-time state-space matrices (x_dot = A x + B u + E d, y = C x)
%   eqns    : struct with human-readable state and output equations (strings)
%   diag    : diagnostics struct (valve slopes, header sensitivities, etc.)
%
% Assumes your init script has populated base workspace:
%   x0_1, x0_2 (states [Vwt;p;a3;Vsd])
%   u0_1, u0_2 (inputs [Q;qf;q4;hf])
%   hdr_IC (struct with p_hdr0, opening1, opening2, mdot1, mdot2, mdot_total)
%   Vd, Vr, Adc, b, tau_s, m_metal, m_riser, k1, k2
%
% Requires on path:
%   boiler_dynamics_block.m
%   Header_NIn1Out_MPa_SC_TYPED.m
%   SeriesValvePipe_HEM_DW_SC.m
%   valve_boiler_to_header.m

  if nargin < 1
    useSeriesValvePipe = true;
  end

  % --- Pull operating points and physical params from base workspace ---
  x0_1   = evalin('base','x0_1');     % [Vwt;p;a3;Vsd]
  x0_2   = evalin('base','x0_2');
  u0_1   = evalin('base','u0_1');     % [Q;qf;q4;hf]
  u0_2   = evalin('base','u0_2');
  hdr_IC = evalin('base','hdr_IC');   % struct: p_hdr0, opening1, opening2, mdot1, mdot2, mdot_total
  Vd     = evalin('base','Vd');

  params = struct('b',evalin('base','b'), 'tau_s',evalin('base','tau_s'), ...
                  'm_metal',evalin('base','m_metal'), 'm_riser',evalin('base','m_riser'), ...
                  'Vr',evalin('base','Vr'), 'Adc',evalin('base','Adc'));
  params1 = params; params1.k = evalin('base','k1');
  params2 = params; params2.k = evalin('base','k2');

  % Header operating point (adjust header quality and enthalpies if you have them)
  hdr_op.p_hdr = hdr_IC.p_hdr0;
  hdr_op.x_hdr = 0.10;                 % header quality (assumed)
  hdr_op.m_out = hdr_IC.mdot_total;    % disturbance base (consumer draw)
  hdr_op.m_in  = [hdr_IC.mdot1; hdr_IC.mdot2];
  hdr_op.h_in  = [u0_1(4); u0_2(4)];   % feed enthalpy placeholder or valve outlet enthalpies

  % --- Linearize boilers and header blocks ---
  level1 = @(x) x(1)/Vd; level2 = @(x) x(1)/Vd;
  [A1,B1,~,dLdx1] = lin_boiler_at_op(x0_1,u0_1,params1,hdr_op.p_hdr,level1);
  [A2,B2,~,dLdx2] = lin_boiler_at_op(x0_2,u0_2,params2,hdr_op.p_hdr,level2);

  [Ah,Eh,Sh_in,p_out_slope,~] = lin_header_at_op(hdr_op, struct('D',0.40,'L',50.0), ...
                                                  struct('rough_m',4.5e-5,'K_minor',0.5,'mu_Pa_s',1.5e-5));

  % --- Valve sensitivities at the operating point (choose your block) ---
  opening1 = hdr_IC.opening1; opening2 = hdr_IC.opening2;
  x_up1 = x0_1(3); x_up2 = x0_2(3); % upstream quality proxies (adjust if needed)

  if useSeriesValvePipe
    s1 = lin_series_valve_pipe_sens(x0_1(2), x_up1, hdr_op.p_hdr, opening1);
    s2 = lin_series_valve_pipe_sens(x0_2(2), x_up2, hdr_op.p_hdr, opening2);
    dmdot_dp_drum_1 = s1.dmdot_dp_drum; dmdot_dp_hdr_1 = s1.dmdot_dp_hdr;
    dmdot_dp_drum_2 = s2.dmdot_dp_drum; dmdot_dp_hdr_2 = s2.dmdot_dp_hdr;
  else
    s1 = lin_valve_boiler_to_header_sens(x0_1(2), hdr_op.p_hdr, opening1, x_up1);
    s2 = lin_valve_boiler_to_header_sens(x0_2(2), hdr_op.p_hdr, opening2, x_up2);
    dmdot_dp_drum_1 = s1.dmdot_dp_up;  dmdot_dp_hdr_1 = s1.dmdot_dp_dn;
    dmdot_dp_drum_2 = s2.dmdot_dp_up;  dmdot_dp_hdr_2 = s2.dmdot_dp_dn;
  end

  % --- Global state ordering: x = [x1(4); x2(4); p_hdr; x_hdr], u = [qf1;qf2;Q1;Q2], d = [m_load] ---
  A = zeros(10,10); B = zeros(10,4); E = zeros(10,1);

  % Place boiler blocks and input maps
  A(1:4,1:4)   = A1;  B(1:4,[1 3]) = B1;   % [qf1,Q1]
  A(5:8,5:8)   = A2;  B(5:8,[2 4]) = B2;   % [qf2,Q2]

  % Header block and disturbance path
  A(9:10,9:10) = Ah;
  E(9:10,1)    = Eh;

  % Header pressure → boiler dynamics via q4 path (chain rule):
  % ∂dxdt/∂p_hdr ≈ (∂dxdt/∂q4)·(∂q4/∂p_hdr)
  dfdq4_1 = numjacobian(@(q4) boiler_dynamics_block(x0_1,[u0_1(1);u0_1(2);q4;u0_1(4)], ...
                                     params1.b,params1.k,params1.Vr,params1.Adc, ...
                                     params1.m_metal,params1.m_riser,params1.tau_s), u0_1(3), 1e-4);
  dfdq4_2 = numjacobian(@(q4) boiler_dynamics_block(x0_2,[u0_2(1);u0_2(2);q4;u0_2(4)], ...
                                     params2.b,params2.k,params2.Vr,params2.Adc, ...
                                     params2.m_metal,params2.m_riser,params2.tau_s), u0_2(3), 1e-4);

  A(1:4,9) = dfdq4_1 * dmdot_dp_hdr_1 * p_out_slope; % small-signal p_out ≈ p_hdr
  A(5:8,9) = dfdq4_2 * dmdot_dp_hdr_2 * p_out_slope;

  % Boiler drum pressures → header states via inflow sensitivities:
  % Sh_in maps ∂[dpdt;dxdt]/∂(m_in_k). Chain with ∂m_in_k/∂p_drum_k:
  A(9:10,2) = Sh_in(:,1) * dmdot_dp_drum_1;  % p1 → header states
  A(9:10,6) = Sh_in(:,2) * dmdot_dp_drum_2;  % p2 → header states

  % Outputs: y = [L1; L2; p_hdr]
  C = zeros(3,10);
  C(1,1:4) = dLdx1;    % level 1 linearization
  C(2,5:8) = dLdx2;    % level 2 linearization
  C(3,9)   = 1;        % header pressure

  % Pretty-print equations
  state_names  = {'Vwt1','p1','a3_1','Vsd1','Vwt2','p2','a3_2','Vsd2','p_hdr','x_hdr'};
  input_names  = {'qf1','qf2','Q1','Q2'};
  dist_names   = {'m_load'};
  output_names = {'L1','L2','p_hdr'};
  eqns = print_state_equations(A,B,E,C,state_names,input_names,dist_names,output_names);

  % Diagnostics for auditing
  diag = struct();
  diag.valve1 = s1;
  diag.valve2 = s2;
  diag.p_out_slope = p_out_slope;
  diag.Sh_in = Sh_in;
  diag.header_op = hdr_op;
end

% ========================================================================
% Boiler linearisation (single unit)
% ========================================================================
function [Ai,Bi,Sip,dLdx] = lin_boiler_at_op(x0,u0,params,p_hdr0,level_fn)
% x0=[Vwt;p;a3;Vsd], u0=[Q;qf;q4;hf]
  f_dyn = @(x,u) boiler_dynamics_block(x,u,params.b,params.k,params.Vr,params.Adc, ...
                                       params.m_metal,params.m_riser,params.tau_s);

  Ai = numjacobian(@(x) f_dyn(x,u0), x0, 1e-5);
  Bqf = numjacobian(@(qf) f_dyn(x0,[u0(1);qf;u0(3);u0(4)]), u0(2), 1e-5);
  BQ  = numjacobian(@(Q ) f_dyn(x0,[Q;u0(2);u0(3);u0(4)]),   u0(1), 1e-5);
  Bi  = [Bqf,BQ];

  % Header-pressure sensitivity will be injected in assembly via valve slopes
  Sip = zeros(4,1);

  % Output level gradient
  dLdx = numjacobian(@(x) level_fn(x), x0, 1e-6).';
end

% ========================================================================
% Header linearisation (dynamic 2-state: pressure, quality)
% ========================================================================
function [Ah,Eh,Sh_in,p_out_slope,Ch] = lin_header_at_op(hdr_op,geom,loss)
  f_hdr = @(p,x,m_out,m_in,h_in) Header_NIn1Out_MPa_SC_TYPED(p,x,m_out,m_in,h_in,geom,loss);

  p0=hdr_op.p_hdr; x0=hdr_op.x_hdr; mout=hdr_op.m_out; minv=hdr_op.m_in; hinv=hdr_op.h_in;

  Jx = numjacobian(@(vx) hdr_pack(f_hdr,vx(1),vx(2),mout,minv,hinv), [p0;x0], 1e-5);
  Ah = Jx;                      % ∂[dpdt;dxdt]/∂[p_hdr;x_hdr]

  Eh = numjacobian(@(m) hdr_pack(f_hdr,p0,x0,m,minv,hinv), mout, 1e-5); % disturbance m_out

  nin=numel(minv); Sh_in=zeros(2,nin);
  for k=1:nin
      mk_fun=@(mk) hdr_pack(f_hdr,p0,x0,mout,replace_at(minv,k,mk),hinv);
      Sh_in(:,k)=numjacobian(mk_fun,minv(k),1e-5);
  end

  p_out_slope = numjacobian(@(ph) hdr_pout(f_hdr,ph,x0,mout,minv,hinv), p0, 1e-5);
  Ch=[1,0];
end

function y=hdr_pack(f_hdr,p,x,m_out,m_in,h_in)
  [dpdt,dxdt]=f_hdr(p,x,m_out,m_in,h_in); y=[dpdt;dxdt];
end
function p_out=hdr_pout(f_hdr,p,x,m_out,m_in,h_in)
  [~,~,~,p_out]=f_hdr(p,x,m_out,m_in,h_in);
end
function v=replace_at(v,k,val)
  v=v(:); v(k)=val;
end

% ========================================================================
% Valve sensitivity: SeriesValvePipe_HEM_DW_SC
% ========================================================================
function sens = lin_series_valve_pipe_sens(p_drum, x_up, p_hdr, opening)
% Numeric sensitivities of mdot to drum/header pressures and opening
% using SeriesValvePipe_HEM_DW_SC.
  base = @(pd,ph,opn) SeriesValvePipe_HEM_DW_SC(pd, x_up, ph, opn);

  % Base evaluation
  [mdot0, ~, dpv0, dpp0, ~, ~, ~] = base(p_drum, p_hdr, opening);

  % Steps
  dp = 1e-4;         % MPa
  do = 1e-3;         % opening fraction

  % Drum pressure slope
  [mdot_p, ~, ~, ~, ~, ~, ~] = base(p_drum + dp, p_hdr, opening);
  [mdot_m, ~, ~, ~, ~, ~, ~] = base(p_drum - dp, p_hdr, opening);
  dmdot_dp_drum = (mdot_p - mdot_m) / (2*dp);

  % Header pressure slope
  [mdot_p, ~, ~, ~, ~, ~, ~] = base(p_drum, p_hdr + dp, opening);
  [mdot_m, ~, ~, ~, ~, ~, ~] = base(p_drum, p_hdr - dp, opening);
  dmdot_dp_hdr = (mdot_p - mdot_m) / (2*dp);

  % Opening slope
  [mdot_p, ~, ~, ~, ~, ~, ~] = base(p_drum, p_hdr, min(max(opening + do,0),1));
  [mdot_m, ~, ~, ~, ~, ~, ~] = base(p_drum, p_hdr, min(max(opening - do,0),1));
  dmdot_dopening = (mdot_p - mdot_m) / (2*do);

  sens = struct('dmdot_dp_drum', dmdot_dp_drum, ...
                'dmdot_dp_hdr',  dmdot_dp_hdr,  ...
                'dmdot_dopening',dmdot_dopening,...
                'mdot_op',       mdot0,         ...
                'dp_valve_op',   dpv0,          ...
                'dp_pipe_op',    dpp0);
end

% ========================================================================
% Valve sensitivity: valve_boiler_to_header
% ========================================================================
function sens = lin_valve_boiler_to_header_sens(p_up, p_dn, opening, x_boiler)
% Numeric sensitivities using valve_boiler_to_header.
  base = @(pu,pd,opn) valve_boiler_to_header(pu, pd, opn, x_boiler);

  % Base eval
  [mdot0, ~, regime0, ~, dpv0] = base(p_up, p_dn, opening);

  % Steps
  dp = 1e-4;  do = 1e-3;

  % Upstream (drum) pressure slope
  [mdot_p, ~, ~, ~, ~] = base(p_up + dp, p_dn, opening);
  [mdot_m, ~, ~, ~, ~] = base(p_up - dp, p_dn, opening);
  dmdot_dp_up = (mdot_p - mdot_m) / (2*dp);

  % Downstream (header) pressure slope
  [mdot_p, ~, ~, ~, ~] = base(p_up, p_dn + dp, opening);
  [mdot_m, ~, ~, ~, ~] = base(p_up, p_dn - dp, opening);
  dmdot_dp_dn = (mdot_p - mdot_m) / (2*dp);

  % Opening slope
  [mdot_p, ~, ~, ~, ~] = base(p_up, p_dn, min(max(opening + do,0),1));
  [mdot_m, ~, ~, ~, ~] = base(p_up, p_dn, min(max(opening - do,0),1));
  dmdot_dopening = (mdot_p - mdot_m) / (2*do);

  sens = struct('dmdot_dp_up',    dmdot_dp_up, ...
                'dmdot_dp_dn',    dmdot_dp_dn, ...
                'dmdot_dopening', dmdot_dopening, ...
                'mdot_op',        mdot0, ...
                'dp_valve_op',    dpv0,  ...
                'regime',         regime0);
end

% ========================================================================
% Numeric Jacobian (central difference)
% ========================================================================
function J = numjacobian(f, x0, relstep)
% Central-difference Jacobian J = ∂f/∂x at x0
  if nargin<3, relstep=1e-5; end
  x0 = x0(:);
  y0 = f(x0);
  m  = numel(y0);
  n  = numel(x0);
  J  = zeros(m,n);
  for k = 1:n
    xk = x0(k);
    delta = max(abs(xk)*relstep, relstep);
    xp = x0; xm = x0;
    xp(k) = xk + delta; xm(k) = xk - delta;
    yp = f(xp); ym = f(xm);
    J(:,k) = (yp - ym) / (2*delta);
  end
end

% ========================================================================
% Equation printer (human-readable state and output equations)
% ========================================================================
function eqns = print_state_equations(A,B,E,C,state_names,input_names,dist_names,output_names)
% Returns a struct of strings and prints them
  nx=size(A,1); nu=size(B,2); nd=size(E,2); ny=size(C,1);
  eqns = struct('states', {cell(nx,1)}, 'outputs', {cell(ny,1)});
  for i=1:nx
    s = sprintf('dx(%s)/dt =', state_names{i});
    for j=1:nx, if abs(A(i,j))>1e-10, s=[s, sprintf(' %+0.5g*%s',A(i,j),state_names{j})]; end, end
    for j=1:nu, if abs(B(i,j))>1e-10, s=[s, sprintf(' %+0.5g*%s',B(i,j),input_names{j})]; end, end
    for j=1:nd, if abs(E(i,j))>1e-10, s=[s, sprintf(' %+0.5g*%s',E(i,j),dist_names{j})]; end, end
    disp(s); eqns.states{i}=s;
  end
  for i=1:ny
    s = sprintf('y(%s) =', output_names{i});
    for j=1:nx, if abs(C(i,j))>1e-10, s=[s, sprintf(' %+0.5g*%s',C(i,j),state_names{j})]; end, end
    disp(s); eqns.outputs{i}=s;
  end
end





%% ========================================================================
%  Test script for two‑boiler + header linearisation
%  Save as test_linearisation.m
%  Run AFTER your init script has populated the base workspace
%  ------------------------------------------------------------------------


disp('>>> Running linearisation test...');

%% 1. Choose valve model
% true  -> use SeriesValvePipe_HEM_DW_SC
% false -> use valve_boiler_to_header
useSeriesValvePipe = true;

%% 2. Call the linearisation function
[A,B,E,C,eqns,diag] = linearise_two_boiler_header_full(useSeriesValvePipe);

%% 3. Inspect matrix sizes
fprintf('\n--- Matrix dimensions ---\n');
fprintf('A: %dx%d\n', size(A));
fprintf('B: %dx%d\n', size(B));
fprintf('E: %dx%d\n', size(E));
fprintf('C: %dx%d\n', size(C));

%% 4. Show the constructed matrices
fprintf('\n--- Matrix A ---\n'); disp(A);
fprintf('\n--- Matrix B ---\n'); disp(B);
fprintf('\n--- Matrix E ---\n'); disp(E);
fprintf('\n--- Matrix C ---\n'); disp(C);

%% 5. Check steady‑state residuals
x0   = [evalin('base','x0_1'); evalin('base','x0_2'); ...
        evalin('base','hdr_IC.p_hdr0'); 0.10];   % header quality guess
u    = [evalin('base','u0_1(2)'); evalin('base','u0_2(2)'); ...
        evalin('base','u0_1(1)'); evalin('base','u0_2(1)')];
d    = evalin('base','hdr_IC.mdot_total');

res  = A*x0 + B*u + E*d;
fprintf('\n--- Steady‑state residual norm ---\n');
disp(norm(res));

%% 6. Print a few equations
fprintf('\n--- Example state equations ---\n');
disp(eqns.states{1});
disp(eqns.states{2});
fprintf('\n--- Example output equations ---\n');
disp(eqns.outputs{1});
disp(eqns.outputs{3});

%% 7. Build state‑space object and simulate
sys = ss(A,B,C,zeros(size(C,1),size(B,2)));

t = 0:1:200;                 % 200s horizon
u_in = zeros(length(t),4);   % 4 inputs: [qf1,qf2,Q1,Q2]
u_in(:,3) = 1;               % step of +1 MW on Q1

[y,t] = lsim(sys,u_in,t);

%% 8. Plot results
figure;
plot(t,y);
xlabel('Time [s]');
ylabel('Outputs');
legend('L1','L2','p_{hdr}');
title('Linearised model response to +1 MW step on Q1');

%% 9. Diagnostics
fprintf('\n--- Valve 1 sensitivities ---\n');
disp(diag.valve1);
fprintf('\n--- Valve 2 sensitivities ---\n');
disp(diag.valve2);

% Build state-space model from your matrices
sys_all = ss(A,B,C,0);

% Extract header pressure output (row 3 of C) and Q1 input (column 3 of B)
sys_hdr   = ss(A,B,C(3,:),0);   % header pressure as output
G_Q1_to_p = sys_hdr(:,3);       % input #3 = Q1