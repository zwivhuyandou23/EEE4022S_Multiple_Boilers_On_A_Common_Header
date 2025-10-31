function [p_hdr0, info] = solve_hdr_for_total_qs_two_boilers( ...
    p_drum1_MPa, x_up1, opening1, ...
    p_drum2_MPa, x_up2, opening2, ...
    qs_total_target)
% Solve header pressure p_hdr0 such that:
%   mdot1(p_hdr) + mdot2(p_hdr) = qs_total_target
% Boiler 1: SeriesValvePipe_HEM_DW_SC (valve + pipe)
% Boiler 2: valve_boiler_to_header (valve only)
%
% Inputs:
%   p_drum1_MPa, x_up1, opening1 : Boiler 1 upstream state
%   p_drum2_MPa, x_up2, opening2 : Boiler 2 upstream state
%   qs_total_target               : desired total inflow [kg/s]
%
% Outputs:
%   p_hdr0 : solved header pressure [MPa]
%   info   : struct with diagnostic results

% --- Clamp inputs ---
p_drum1_MPa = max(0.1, double(p_drum1_MPa));
p_drum2_MPa = max(0.1, double(p_drum2_MPa));
x_up1       = max(0.0, min(1.0, double(x_up1)));
x_up2       = max(0.0, min(1.0, double(x_up2)));
opening1    = max(0.0, min(1.0, double(opening1)));
opening2    = max(0.0, min(1.0, double(opening2)));
qs_total_target = max(1e-6, double(qs_total_target));

% --- Residual function: total_flow(p_hdr) - target ---
function r = resid(p_hdr)
    [mdot1, ~, ~, ~] = SeriesValvePipe_HEM_DW_SC(p_drum1_MPa, x_up1, p_hdr, opening1);
    [mdot2, ~, ~, ~, ~, ~] = valve_boiler_to_header(p_drum2_MPa, p_hdr, opening2, x_up2);
    r = (mdot1 + mdot2) - qs_total_target;
end

% --- Bracket search for feasible range ---
p_min = 0.10;  % MPa
p_max = min(p_drum1_MPa, p_drum2_MPa) - 1e-3;
grid  = linspace(p_min, p_max, 60);
vals  = arrayfun(@resid, grid);
flows = vals + qs_total_target;

info = struct();
info.p_grid = grid;
info.total_flow_grid = flows;

% --- Feasibility check ---
if max(flows) < qs_total_target - 1e-6
    % Even at lowest header pressure, cannot reach target
    p_hdr0 = p_min;
    info.feasible = false;
    info.max_flow = max(flows);
    return;
end

% --- Find crossing interval ---
idx = find(vals(1:end-1).*vals(2:end) <= 0, 1, 'first');
if isempty(idx)
    % If no exact crossing, pick nearest
    [~,iBest] = min(abs(vals));
    p_hdr0 = grid(iBest);
else
    p_lower = grid(idx);
    p_upper = grid(idx+1);
    p_hdr0  = fzero(@resid, [p_lower, p_upper]);
end

% --- Safety clamp ---
p_hdr0 = max(p_min, min(p_hdr0, p_max));

% --- Diagnostics at solution ---
[mdot1, h_out1, dp_valve1, dp_pipe1] = ...
    SeriesValvePipe_HEM_DW_SC(p_drum1_MPa, x_up1, p_hdr0, opening1);

[mdot2, h_out2, regime2, Tsat2, dp_valve2, x2] = ...
    valve_boiler_to_header(p_drum2_MPa, p_hdr0, opening2, x_up2);

info.feasible     = true;
info.mdot1        = mdot1;
info.mdot2        = mdot2;
info.total_flow   = mdot1 + mdot2;
info.h_out1       = h_out1;
info.h_out2       = h_out2;
info.dp_valve1    = dp_valve1;
info.dp_pipe1     = dp_pipe1;
info.dp_valve2    = dp_valve2;
info.regime2      = regime2;
info.Tsat2        = Tsat2;
info.x2           = x2;
end