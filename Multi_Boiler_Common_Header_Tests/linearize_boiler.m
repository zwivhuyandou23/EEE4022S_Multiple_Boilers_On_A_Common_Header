function [A,B,C,D,sys,G] = linearize_boiler(x0,u0,b,k,Vr,Adc,m_metal,m_riser,tau_s,Vdc,Ad)
% LINEARIZE_BOILER_SS
%   Linearises the 4‑state drum‑boiler model around operating point (x0,u0).
%   States: x = [Vwt; p_MPa; a3; Vsd]
%   Inputs: u = [Q_MW; qf; q4; hf]
%   Outputs: y = [drum pressure (MPa); drum water level (m)]
%
%   Returns:
%     A,B,C,D  - state-space matrices
%     sys      - state-space object
%     G        - transfer function object (MIMO)

    n_states = length(x0);
    n_inputs = length(u0);

    A = zeros(n_states,n_states);
    B = zeros(n_states,n_inputs);

    % --- Nominal derivative at operating point
    [f0,~] = boiler_dynamics_block(x0,u0,b,k,Vr,Adc,m_metal,m_riser,tau_s);

    % --- A matrix (df/dx) ---
    for i = 1:n_states
        dx = 1e-6*max(abs(x0(i)),1);   % scaled perturbation
        x_plus = x0; x_plus(i) = x_plus(i) + dx;
        x_minus= x0; x_minus(i)= x_minus(i) - dx;

        f_plus = boiler_dynamics_block(x_plus,u0,b,k,Vr,Adc,m_metal,m_riser,tau_s);
        f_minus= boiler_dynamics_block(x_minus,u0,b,k,Vr,Adc,m_metal,m_riser,tau_s);

        A(:,i) = (f_plus - f_minus)/(2*dx);
    end

    % --- B matrix (df/du) ---
    for j = 1:n_inputs
        du = 1e-6*max(abs(u0(j)),1);   % scaled perturbation
        u_plus = u0; u_plus(j) = u_plus(j) + du;
        u_minus= u0; u_minus(j)= u0(j) - du;

        f_plus = boiler_dynamics_block(x0,u_plus,b,k,Vr,Adc,m_metal,m_riser,tau_s);
        f_minus= boiler_dynamics_block(x0,u_minus,b,k,Vr,Adc,m_metal,m_riser,tau_s);

        B(:,j) = (f_plus - f_minus)/(2*du);
    end

    % --- Output matrices ---
    % y1 = drum pressure = p_MPa (state 2)
    C1 = [0 1 0 0];

    % y2 = drum water level = (Vwt - Vdc - (1-av)*Vr + Vsd)/Ad
    % Need void fraction sensitivities at op point
    [~,~,dav_dp,dav_da3] = get_boiler_calcs(x0(2),x0(3));

    C2 = [ 1/Ad , (Vr/Ad)*dav_dp , (Vr/Ad)*dav_da3 , 1/Ad ];

    C = [C1; C2];
    D = zeros(2,n_inputs);

    % --- Build state-space and transfer function objects ---
    sys = ss(A,B,C,D);
    G   = tf(sys);   % MIMO transfer function (2 outputs × 4 inputs)

end