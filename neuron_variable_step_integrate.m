function [vMembrane, n, m, h, iLeak, iK, iNa, timeStamps] = neuron_variable_step_integrate(...
    iInput, restPotential, eqPotLeak, eqPotK, eqPotNa, gLeak, gKMax, gNaMax, ...
    capacitance, simulationTime, timeStepSize)
% Description: ODE15s integration of the Hodgkin-Huxley model.
% Units: Time (ms), Current (uA/cm^2), Potential (mV), Conductance (mS/cm^2)
% Note: timeStepSize defines the output evaluation points, not the solver step size.

%% 1. Simulation Setup
% totalNumSteps = simulationTime / timeStepSize;
% dt = timeStepSize;

%% 2. Gating Variable Dynamics (Hodgkin-Huxley)
% V is relative to resting membrane potential (this is the convention in the hodgkin huxley paper)
V = @(vMem) -(vMem - restPotential);

% Potassium (n)
get_n_alpha = @(vMem) 0.01 * (10 + V(vMem)) / (exp((10 + V(vMem)) / 10) - 1);
get_n_beta  = @(vMem) 0.125 * exp(V(vMem) / 80);
get_n_tau   = @(vMem) 1 / (get_n_alpha(vMem) + get_n_beta(vMem));
get_n_ss    = @(vMem) get_n_alpha(vMem) / (get_n_alpha(vMem) + get_n_beta(vMem));

% Sodium Activation (m)
get_m_alpha = @(vMem) 0.1 * (25 + V(vMem)) / (exp((25 + V(vMem)) / 10) - 1);
get_m_beta  = @(vMem) 4 * exp(V(vMem) / 18);
get_m_tau   = @(vMem) 1 / (get_m_alpha(vMem) + get_m_beta(vMem));
get_m_ss    = @(vMem) get_m_alpha(vMem) / (get_m_alpha(vMem) + get_m_beta(vMem));

% Sodium Inactivation (h)
get_h_alpha = @(vMem) 0.07 * exp(V(vMem) / 20);
get_h_beta  = @(vMem) 1 / (exp((30 + V(vMem)) / 10) + 1);
get_h_tau   = @(vMem) 1 / (get_h_alpha(vMem) + get_h_beta(vMem));
get_h_ss    = @(vMem) get_h_alpha(vMem) / (get_h_alpha(vMem) + get_h_beta(vMem));

%% Freaky matlab ode function declaration:
    function dfdt = odesys(t,state)
        vMembrane_ = state(1);
        n_ = state(2);
        m_ = state(3);
        h_ = state(4);

        ileak_ = gLeak  * (vMembrane_ - eqPotLeak);
        iK_  = gKMax  *      n_^4 * (vMembrane_ - eqPotK);
        iNa_ = gNaMax * m_^3 * h_ * (vMembrane_ - eqPotNa);

        dV = (1 / capacitance) * (iInput(t) - ileak_ - iK_ - iNa_);
        dn = (get_n_ss(vMembrane_) - n_) / get_n_tau(vMembrane_);
        dm = (get_m_ss(vMembrane_) - m_) / get_m_tau(vMembrane_);
        dh = (get_h_ss(vMembrane_) - h_) / get_h_tau(vMembrane_);

        dfdt = [dV; dn; dm; dh];

    end

% Initial Conditions
% V, nprob, mprob, hprob
state0 = [restPotential, get_n_ss(restPotential), get_m_ss(restPotential), get_h_ss(restPotential)];

%% Solve
tSpan = 0: timeStepSize: simulationTime;
options = odeset('RelTol',1e-6, 'AbsTol', 1e-8);
% Perform the numerical integration using ODE solver
[timeStamps, solution] = ode15s(@odesys, tSpan, state0, options);

    % Fix output dims
    timeStamps = timeStamps';
    vMembrane  = solution(:, 1)';
    n          = solution(:, 2)';
    m          = solution(:, 3)';
    h          = solution(:, 4)';

        %% Calculate Currents for compatability
    iLeak = gLeak  *           n.^4 .* (vMembrane - eqPotLeak);
    iK    = gKMax  *           n.^4 .* (vMembrane - eqPotK);
    iNa   = gNaMax * m.^3 .* h     .* (vMembrane - eqPotNa);
end