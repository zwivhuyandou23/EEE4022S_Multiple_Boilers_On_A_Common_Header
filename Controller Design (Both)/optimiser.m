function [p_sp, w, Qsp, Loss, Loss_i] = optimiser( ...
     h_in, h_out, m_dot, Q_star, dQ, prev_w)
% OPTIMISER  Supervisory optimiser (2 boilers) using provided enthalpies.
% Loss model is strictly Loss_i = m_dot(i) * (h_in(i) - h_out).
% No enthalpy calculations inside; h_in and h_out are inputs.
%
% Inputs:
%   p_hdr  [MPa]   measured header pressure (used only to keep p_sp in band)
%   h_in   [2x1]   inlet enthalpy per boiler (kJ/kg), provided externally
%   h_out  [1x1]   header outlet enthalpy (kJ/kg), provided externally
%   m_dot  [2x1]   mass flow per boiler to the header (kg/s)
%   Q_star [2x1]   nominal firing points per boiler (MW)
%   dQ     [1x1]   bias from master PI (MW) to be allocated
%   prev_w [2x1]   previous allocation weights (for initialisation)
%
% Outputs:
%   p_sp   [MPa]   header pressure setpoint (clamped within band)
%   w      [2x1]   allocation weights (sum = 1)
%   Qsp    [2x1]   firing setpoints: Qsp = Q_star + w*dQ (MW)
%   Loss   [kW]    total power loss in header (kJ/s ≈ kW)
%   Loss_i [2x1]   per‑boiler power loss (kJ/s ≈ kW)

% --- Self-contained constants (tunable band and bounds) ---
p_nom = 8.5;           % [MPa] nominal header pressure setpoint
p_min = 8.501;           % [MPa] minimum allowed header pressure
p_max = 8.499;           % [MPa] maximum allowed header pressure
w_min = [0.0; 0.0];    % minimum share for [Boiler 1; Boiler 2]
w_max = [1.0; 1.0];    % maximum share for [Boiler 1; Boiler 2]

% --- Strict loss model (no extra penalties) ---
Loss = m_dot .* (h_in - h_out);  % kJ/s (≈ kW)


% --- Choose allocation weights to bias toward lower-loss branch ---
% Heuristic objective: J = w' * Loss_i (minimise weighted loss)
w1_grid = linspace(w_min(1), w_max(1), 101);
bestJ = inf; bestw = prev_w;
for w1 = w1_grid
    w2 = 1 - w1;
    if w2 < w_min(2) || w2 > w_max(2), continue; end
    w_try = [w1; w2];
    J = w_try.' * Loss_i;  % lower J favors assigning ΔQ to lower-loss branch
    if J < bestJ
        bestJ = J;
        bestw = w_try;
    end
end
w = bestw;

% --- Dispatch firing setpoints (allocator relationship) ---


% --- Header pressure setpoint (kept within band) ---
p_sp = max(p_min, min(p_nom, p_max));

end

%% === Loss-only optimiser demo (bias allocation toward lower-loss branch) ===
% Build loss inputs from init (h_in provided by boiler models; here use saturated vapor enthalpy)
pbar1 = x0_1(2) * 10;  % Boiler 1 (far) drum pressure [bar]
pbar2 = x0_2(2) * 10;  % Boiler 2 (close) drum pressure [bar]
h_s1  = XSteam('hV_p', pbar1);
h_s2  = XSteam('hV_p', pbar2);
h_out = XSteam('hV_p', p_hdr_star*10);   % header enthalpy (dry saturated)

h_in  = [h_s1; h_s2];         % [far; close]
m_dot = [u0_1(3); u0_2(3)];   % q4 [far; close]
Q_star= [u0_1(1); u0_2(1)];   % Q [far; close]
dQ    = 0.0;                  % demo (PI output at this instant)
prev_w= [w1; w2];

[p_sp, w, Qsp, Loss, Loss_i] = optimiser_loss_only(p_hdr_star, h_in, h_out, m_dot, Q_star, dQ, prev_w);

fprintf('\n=== Optimiser (loss-only) ===\n');
fprintf('p_sp=%.3f MPa, w=[%.3f %.3f], Qsp=[%.3f %.3f] MW\n', p_sp, w(1), w(2), Qsp(1), Qsp(2));
fprintf('Loss total=%.2f kW, Loss_i=[%.2f %.2f] kW\n', Loss, Loss_i(1), Loss_i(2));
