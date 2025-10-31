function props = get_steam_props(p_bar)
% This is the SINGLE SOURCE OF TRUTH for steam property calculations.
% It is an external M-file.

coder.extrinsic('XSteam');

% Pre-initialize the structure for the Coder
props=struct('rhoL',0.0,'rhoV',0.0,'hL',0.0,'hV',0.0,'vL',0.0,'vV',0.0,...
             'drhoL_dpbar',0.0,'drhoV_dpbar',0.0,'dhL_dpbar',0.0,'dhV_dpbar',0.0,...
             'dvL_dpbar',0.0,'dvV_dpbar',0.0,'Tsat',0.0,'dTs_dpbar',0.0);

% Call XSteam and cast outputs to double
props.rhoL=double(XSteam('rhoL_p',p_bar));
props.rhoV=double(XSteam('rhoV_p',p_bar));
props.hL=double(XSteam('hL_p',p_bar));
props.hV=double(XSteam('hV_p',p_bar));
props.vL=double(XSteam('vL_p',p_bar));
props.vV=double(XSteam('vV_p',p_bar));
props.Tsat=double(XSteam('Tsat_p',p_bar));

% Calculate derivatives
dp_step=1e-5;
p_bar_plus=p_bar+dp_step;
p_bar_minus=p_bar-dp_step;

props.drhoL_dpbar=(double(XSteam('rhoL_p',p_bar_plus))-double(XSteam('rhoL_p',p_bar_minus)))/(2*dp_step);
props.drhoV_dpbar=(double(XSteam('rhoV_p',p_bar_plus))-double(XSteam('rhoV_p',p_bar_minus)))/(2*dp_step);
props.dhL_dpbar=(double(XSteam('hL_p',p_bar_plus))-double(XSteam('hL_p',p_bar_minus)))/(2*dp_step);
props.dhV_dpbar=(double(XSteam('hV_p',p_bar_plus))-double(XSteam('hV_p',p_bar_minus)))/(2*dp_step);
props.dvL_dpbar=(double(XSteam('vL_p',p_bar_plus))-double(XSteam('vL_p',p_bar_minus)))/(2*dp_step);
props.dvV_dpbar=(double(XSteam('vV_p',p_bar_plus))-double(XSteam('vV_p',p_bar_minus)))/(2*dp_step);
props.dTs_dpbar=(double(XSteam('Tsat_p',p_bar_plus))-double(XSteam('Tsat_p',p_bar_minus)))/(2*dp_step);
end