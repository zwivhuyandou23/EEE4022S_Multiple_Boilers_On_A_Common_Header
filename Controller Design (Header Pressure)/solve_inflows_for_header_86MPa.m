function result = solve_inflows_for_header_86MPa( ...
        p_hdr_MPa, ...
        p1_MPa, x1_up, ...
        p2_MPa, x2_up, ...
        m_out_kg_s, ...
        w1, w2)
% solve_inflows_for_header_86MPa
% Compute mass inflows from two boilers such that the header is held at
% p_hdr_MPa (e.g. 8.6 MPa) in steady state, assuming equal inlet enthalpies
% to the header. Returns valve openings and inflows that satisfy:
%   mdot1 + mdot2 = m_out_kg_s
% Optionally enforces a flow split via weights w1:w2 (default: design ratio).
%
% Inputs:
%   p_hdr_MPa   target header pressure [MPa], e.g., 8.6
%   p1_MPa      boiler 1 drum pressure [MPa]
%   x1_up       boiler 1 upstream quality [-]
%   p2_MPa      boiler 2 drum pressure [MPa]
%   x2_up       boiler 2 upstream quality [-]
%   m_out_kg_s  downstream outflow from header [kg/s]
%   w1, w2      nonnegative weights for flow split (optional)
%
% Output struct fields:
%   opening1, opening2   valve openings [0..1]
%   mdot1, mdot2         inflows [kg/s]
%   mdot_total           sum inflow [kg/s] (≈ m_out_kg_s)
%   p_hdr_MPa            header pressure setpoint (echo)
%   share                mdot1/(mdot1+mdot2)

% Defaults: keep original design split if not provided
if nargin < 7 || isempty(m_out_kg_s), error('m_out_kg_s required'); end
if nargin < 8 || isempty(w1), w1 = 58.5; end
if nargin < 9 || isempty(w2), w2 = 54.0; end

% Normalize weights to target shares
w1 = max(0, double(w1)); w2 = max(0, double(w2));
if (w1 + w2) <= 0, w1 = 1; w2 = 1; end
share1 = w1 / (w1 + w2);
share2 = 1 - share1;

% Target inflows respecting mass balance and desired split
m1_target = share1 * m_out_kg_s;
m2_target = share2 * m_out_kg_s;

% Monotone opening search via bisection
opening1 = solve_opening_for_flow(@SeriesValvePipe_HEM_DW_SC, ...
    p1_MPa, x1_up, p_hdr_MPa, m1_target);
opening2 = solve_opening_for_flow(@SeriesValvePipe_HEM_DW_SC, ...
    p2_MPa, x2_up, p_hdr_MPa, m2_target);

% Evaluate resulting inflows
[mdot1, ~, ~, ~] = eval_flow(@SeriesValvePipe_HEM_DW_SC, p1_MPa, x1_up, p_hdr_MPa, opening1);
[mdot2, ~, ~, ~] = eval_flow(@SeriesValvePipe_HEM_DW_SC, p2_MPa, x2_up, p_hdr_MPa, opening2);
mdot_total = mdot1 + mdot2;

% If the sum deviates due to nonlinearity, do a small correction loop on openings
tol = 0.25; max_iter = 5;
for k = 1:max_iter
    err = mdot_total - m_out_kg_s;
    if abs(err) < tol, break; end
    % Distribute correction proportional to shares
    m1_target = max(0, m1_target - share1 * err);
    m2_target = max(0, m2_target - share2 * err);
    opening1 = solve_opening_for_flow(@SeriesValvePipe_HEM_DW_SC, p1_MPa, x1_up, p_hdr_MPa, m1_target);
    opening2 = solve_opening_for_flow(@SeriesValvePipe_HEM_DW_SC, p2_MPa, x2_up, p_hdr_MPa, m2_target);
    [mdot1, ~, ~, ~] = eval_flow(@SeriesValvePipe_HEM_DW_SC, p1_MPa, x1_up, p_hdr_MPa, opening1);
    [mdot2, ~, ~, ~] = eval_flow(@SeriesValvePipe_HEM_DW_SC, p2_MPa, x2_up, p_hdr_MPa, opening2);
    mdot_total = mdot1 + mdot2;
end

result = struct();
result.opening1   = opening1;
result.opening2   = opening2;
result.mdot1      = mdot1;
result.mdot2      = mdot2;
result.mdot_total = mdot_total;
result.p_hdr_MPa  = p_hdr_MPa;
result.share      = mdot1 / max(mdot_total, eps);

end

% === Opening solver for a single branch (bisection on [0,1]) ===
function opening = solve_opening_for_flow(fun, p_drum, x_up, p_hdr, m_target)
lo = 0.0; hi = 1.0;
for i = 1:24
    mid = 0.5*(lo + hi);
    m_mid = eval_flow(fun, p_drum, x_up, p_hdr, mid);
    if m_mid < m_target
        lo = mid;
    else
        hi = mid;
    end
end
opening = max(0.0, min(1.0, 0.5*(lo + hi)));
end

% === Evaluate flow for a given opening (wrapper) ===
function [mdot, h_out, dpv, dpp] = eval_flow(fun, p_drum, x_up, p_hdr, opening)
[mdot, h_out, dpv, dpp, ~, ~] = fun(p_drum, x_up, p_hdr, opening);
end