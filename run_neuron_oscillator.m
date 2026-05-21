% ---- Corrected call (all SI: V, S, F, seconds) ----
[L_v, R_v, t] = neuron_oscillator_variable_step_integrate( ...
    -0.060,           ...  % restPotential (V)
    -0.060, 8e-9,     ...  % eqPotLeak (V), gLeak (S)
     0.045, 200e-9, 7e-9, ...  % eqPotNa, gNaFast, gNaPers
    -0.070, 100e-9, 80e-9, 80e-9, ...  % eqPotK, gKDelRect, gKPers, gKFast
     0.135, 5e-9, 3.2e-9, ...  % eqPotCa, gCaFast, gCaSlow
    -0.021, 4e-9,     ...  % eqPotHypPol, gHypPol
    -0.0625, 60e-9,   ...  % eqPotSyn, gSynS
     0.002, 0.011, -0.020,   ...  % tau_rise, tau_decay, h0Thresh (s)  <-- was ms
     5e-10,           ...  % capacitance (F)
     8);                  % simulationTime (s)  <-- was ms; need ~75 s for ~7 cycles

% ---- Plot (output is in V, convert to mV) ----
figure;
subplot(2,1,1);
plot(t, L_v * 1000, 'b', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
title('Left Neuron'); grid on;

subplot(2,1,2);
plot(t, R_v * 1000, 'r', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
title('Right Neuron'); grid on;