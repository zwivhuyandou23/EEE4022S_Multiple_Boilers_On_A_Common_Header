function [props, av, dav_dp, dav_da3] = get_boiler_calcs(p_MPa, a3)
% Computes steam properties and their derivatives at a given pressure and steam quality

coder.extrinsic('XSteam');

% Convert pressure to bar
p_bar = max(1, p_MPa * 10);
dp = 1e-3;  % Small delta for numerical differentiation
a3c = min(max(a3, 1e-6), 1 - 1e-6);  % Clamp a3 to avoid singularities

% Initialize property structure
props = struct('rhoL',0,'rhoV',0,'hL',0,'hV',0,'vL',0,'vV',0,'Tsat',0, ...
               'drhoL_dpbar',0,'drhoV_dpbar',0,'dhL_dpbar',0,'dhV_dpbar',0, ...
               'dvL_dpbar',0,'dvV_dpbar',0,'dTs_dpbar',0);

% Get base properties from XSteam
props.rhoL = double(XSteam('rhoL_p', p_bar));
props.rhoV = double(XSteam('rhoV_p', p_bar));
props.hL   = double(XSteam('hL_p', p_bar));
props.hV   = double(XSteam('hV_p', p_bar));
props.vL   = double(XSteam('vL_p', p_bar));
props.vV   = double(XSteam('vV_p', p_bar));
props.Tsat = double(XSteam('Tsat_p', p_bar));

% Central difference derivatives
props.drhoL_dpbar = (double(XSteam('rhoL_p', p_bar + dp)) - double(XSteam('rhoL_p', p_bar - dp))) / (2 * dp);
props.drhoV_dpbar = (double(XSteam('rhoV_p', p_bar + dp)) - double(XSteam('rhoV_p', p_bar - dp))) / (2 * dp);
props.dhL_dpbar   = (double(XSteam('hL_p', p_bar + dp)) - double(XSteam('hL_p', p_bar - dp))) / (2 * dp);
props.dhV_dpbar   = (double(XSteam('hV_p', p_bar + dp)) - double(XSteam('hV_p', p_bar - dp))) / (2 * dp);
props.dvL_dpbar   = (double(XSteam('vL_p', p_bar + dp)) - double(XSteam('vL_p', p_bar - dp))) / (2 * dp);
props.dvV_dpbar   = (double(XSteam('vV_p', p_bar + dp)) - double(XSteam('vV_p', p_bar - dp))) / (2 * dp);
props.dTs_dpbar   = (double(XSteam('Tsat_p', p_bar + dp)) - double(XSteam('Tsat_p', p_bar - dp))) / (2 * dp);

% Compute average steam volume fraction
av = avbar_corrected(a3c, props.vL, props.vV);

% Derivative of av with respect to a3
da = 1e-5;
av_p = avbar_corrected(min(a3c + da, 1 - 1e-9), props.vL, props.vV);
av_m = avbar_corrected(max(a3c - da, 1e-9), props.vL, props.vV);
dav_da3 = (av_p - av_m) / (2 * da);

% Derivative of av with respect to pressure
vL_p = double(XSteam('vL_p', p_bar + dp)); vV_p = double(XSteam('vV_p', p_bar + dp));
vL_m = double(XSteam('vL_p', p_bar - dp)); vV_m = double(XSteam('vV_p', p_bar - dp));
av_p = avbar_corrected(a3c, vL_p, vV_p);
av_m = avbar_corrected(a3c, vL_m, vV_m);
dav_dp = (av_p - av_m) / (2 * dp);  % per bar
end

function av = avbar_corrected(a3, v_w, v_s)
% Computes average specific volume of steam-water mixture using Eq. (12)
den = v_s * (1 - a3) + v_w * a3;
if den <= 0
    av = 0;
    return;
end
gam = a3 * (v_s - v_w) / den;
if abs(gam) < 1e-9
    logt = 1 - gam / 2 + gam^2 / 3;
else
    logt = log1p(gam) / gam;
end
av = (v_w / den) * logt;
end
