function [vMembrane, n, m, h, iLeak, iK, iNa, timeStamps] = neuron_manual_integrate(...
    iInput, restPotential, eqPotLeak, eqPotK, eqPotNa, gLeak, gKMax, gNaMax, ...
    capacitance, simulationTime, timeStepSize)
% Description: Manual integration of the Hodgkin-Huxley model.
% Units: Time (ms), Current (uA/cm^2), Potential (mV), Conductance (mS/cm^2)

%% 1. Simulation Setup
totalNumSteps = simulationTime / timeStepSize;
dt = timeStepSize;

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

%% 3. Preallocation & Initial Conditions
vMembrane  = zeros(1, totalNumSteps);
n          = zeros(1, totalNumSteps);
m          = zeros(1, totalNumSteps);
h          = zeros(1, totalNumSteps);
iLeak      = zeros(1, totalNumSteps);
iK         = zeros(1, totalNumSteps);
iNa        = zeros(1, totalNumSteps);
timeStamps = zeros(1, totalNumSteps);

% Initial States
vMembrane(1)  = restPotential;
timeStamps(1) = 0.0;
n(1)          = get_n_ss(restPotential);
m(1)          = get_m_ss(restPotential);
h(1)          = get_h_ss(restPotential);

% Initial Currents
iLeak(1) = gLeak * (vMembrane(1) - eqPotLeak);
iK(1)    = gKMax * n(1)^4 * (vMembrane(1) - eqPotK);
iNa(1)   = gNaMax * m(1)^3 * h(1) * (vMembrane(1) - eqPotNa);

%% 4. Numerical Integration (Euler Method)
for j = 2:totalNumSteps
    % Calculate Change in Voltage
    dVMembrane = (1 / capacitance) * ...
        (iInput(timeStamps(j-1)) - iLeak(j-1) - iK(j-1) - iNa(j-1));

    % Update Voltage
    vMembrane(j) = vMembrane(j-1) + (dVMembrane * dt);

    % Update Gating Variables (based on previous voltage)
    n(j) = n(j-1) + (get_n_ss(vMembrane(j-1)) - n(j-1)) / get_n_tau(vMembrane(j-1)) * dt;
    m(j) = m(j-1) + (get_m_ss(vMembrane(j-1)) - m(j-1)) / get_m_tau(vMembrane(j-1)) * dt;
    h(j) = h(j-1) + (get_h_ss(vMembrane(j-1)) - h(j-1)) / get_h_tau(vMembrane(j-1)) * dt;

    % Update Currents (based on current voltage and new gating states)
    iLeak(j) = gLeak * (vMembrane(j) - eqPotLeak);
    iK(j)    = gKMax * n(j)^4 * (vMembrane(j) - eqPotK);
    iNa(j)   = gNaMax * m(j)^3 * h(j) * (vMembrane(j) - eqPotNa);

    % Update Time
    timeStamps(j) = timeStamps(j-1) + dt;
end

end