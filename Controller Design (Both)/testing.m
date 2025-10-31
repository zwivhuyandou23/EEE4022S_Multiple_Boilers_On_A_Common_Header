function [opening1, opening2, info] = ...
    solve_openings_for_targets(p_hdr_MPa, ...
                               p_drum1_MPa, x_up1, mdot1_target, ...
                               p_drum2_MPa, x_up2, mdot2_target)
%SOLVE_OPENINGS_FOR_TARGETS  Find valve openings that achieve target outflows
%
%   Given a fixed header pressure, this function adjusts the valve openings
%   for two boilers so that their outflows match specified targets.
%
%   Boiler 1: flows through SeriesValvePipe_HEM_DW_SC (valve + pipe)
%   Boiler 2: flows through valve_boiler_to_header (valve only)
%
%   Inputs:
%       p_hdr_MPa     - common header pressure [MPa]
%       p_drum1_MPa   - Boiler 1 drum pressure [MPa]
%       x_up1         - Boiler 1 steam quality [-]
%       mdot1_target  - desired Boiler 1 outflow [kg/s]
%       p_drum2_MPa   - Boiler 2 drum pressure [MPa]
%       x_up2         - Boiler 2 steam quality [-]
%       mdot2_target  - desired Boiler 2 outflow [kg/s]
%
%   Outputs:
%       opening1      - solved valve opening fraction for Boiler 1 [0..1]
%       opening2      - solved valve opening fraction for Boiler 2 [0..1]
%       info          - struct with diagnostic results (flows, enthalpies, Δp, etc.)

    % --- Objective function: squared error between actual and target flows ---
    function J = obj(vars)
        % Clamp openings to [0,1]
        o1 = max(0,min(1,vars(1)));
        o2 = max(0,min(1,vars(2)));

        % Compute flows for given openings
        [mdot1, ~, ~, ~] = SeriesValvePipe_HEM_DW_SC(p_drum1_MPa, x_up1, p_hdr_MPa, o1);
        [mdot2, ~, ~, ~, ~, ~] = valve_boiler_to_header(p_drum2_MPa, p_hdr_MPa, o2, x_up2);

        % Objective: sum of squared flow errors
        J = (mdot1 - mdot1_target)^2 + (mdot2 - mdot2_target)^2;
    end

    % --- Initial guess: both valves fully open ---
    vars0 = [1.0, 1.0];

    % --- Optimize openings using Nelder–Mead search (fminsearch) ---
    vars_sol = fminsearch(@obj, vars0);

    % Clamp solutions to [0,1]
    opening1 = max(0,min(1,vars_sol(1)));
    opening2 = max(0,min(1,vars_sol(2)));

    % --- Diagnostics at solved openings ---
    [mdot1, h1, dpv1, dpp1] = SeriesValvePipe_HEM_DW_SC(p_drum1_MPa, x_up1, p_hdr_MPa, opening1);
    [mdot2, h2, regime2, Tsat2, dpv2, x2] = valve_boiler_to_header(p_drum2_MPa, p_hdr_MPa, opening2, x_up2);

    % Collect results in struct
    info = struct();
    info.mdot1     = mdot1;
    info.h1        = h1;
    info.dp_valve1 = dpv1;
    info.dp_pipe1  = dpp1;

    info.mdot2     = mdot2;
    info.h2        = h2;
    info.dp_valve2 = dpv2;
    info.regime2   = regime2;
    info.Tsat2     = Tsat2;
    info.x2        = x2;
end


%% ---------------- Example usage ----------------
p_hdr = 8.0;  % header pressure [MPa]

% Boiler 1 targets
p_drum1 = 8.8; 
x_up1   = 0.13; 
mdot1_target = 58.30;
% Boiler 2 targets
p_drum2 = 8.7; 
x_up2   = 0.13; 
mdot2_target = 53;

% Solve for openings
[opening1, opening2, info] = solve_openings_for_targets( ...
    p_hdr, ...
    p_drum1, x_up1, mdot1_target, ...
    p_drum2, x_up2, mdot2_target);

% Display results
fprintf('Solved openings: Boiler1 = %.3f, Boiler2 = %.3f\n', opening1, opening2);
fprintf('Boiler1 flow = %.2f (target %.2f)\n', info.mdot1, mdot1_target);
fprintf('Boiler2 flow = %.2f (target %.2f)\n', info.mdot2, mdot2_target);