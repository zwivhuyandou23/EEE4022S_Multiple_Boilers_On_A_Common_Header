function p_hdr0 = solve_hdr_for_qs_selfcontained(p_drum_MPa, qs, opening0, x_up0)
% Invert SeriesValvePipe_HEM_DW_SC to find p_hdr0 so mdot == qs.
% Returns the lowest header pressure that can accommodate the desired flow.
%
% Inputs:
%   p_drum_MPa : upstream drum pressure [MPa]
%   qs         : desired mass flow [kg/s]
%   opening0   : valve opening fraction [0..1]
%   x_up0      : upstream quality [-] (use ~1.0 if saturated vapor)
% Output:
%   p_hdr0     : header pressure [MPa]

    % Clamp inputs
    p_drum_MPa = max(0.1, double(p_drum_MPa));
    qs         = max(1e-6, double(qs));
    opening0   = max(0.0, min(1.0, double(opening0)));
    x_up0      = max(0.0, min(1.0, double(x_up0)));

    % Residual function: mdot(p_hdr) - qs
    function r = resid(p_hdr)
        [mdot, ~, ~, ~] = SeriesValvePipe_HEM_DW_SC(p_drum_MPa, x_up0, p_hdr, opening0);
        r = mdot - qs;
    end

    % Scan the feasible range
    p_min = 0.1;                          % MPa
    p_max = p_drum_MPa - 1e-1;            % just below drum pressure
    grid  = linspace(p_min, p_max, 50);
    vals  = arrayfun(@resid, grid);

    % Check feasibility
    if all(vals < 0)
        % Even at lowest header pressure, cannot reach qs
        warning('Infeasible: max flow %.2f < requested %.2f', max(vals+qs), qs);
        p_hdr0 = p_min;  % return lowest possible header pressure
        return;
    end

    % Find interval where residual crosses zero
    idx = find(vals(1:end-1).*vals(2:end) <= 0, 1);
    if isempty(idx)
        % Should not happen if feasible
        p_hdr0 = p_min;
        return;
    end

    % Refine with fzero
    p_lower = grid(idx);
    p_upper = grid(idx+1);
    p_hdr0  = fzero(@resid, [p_lower, p_upper]);

    % Safety clamp
    p_hdr0 = max(p_min, min(p_hdr0, p_max));
end