clc;
clear;

%% Neuron Integration Runner
MEMBRANE_CAPACITANCE = 1;
REST_POTENTIAL_mV    = -65;
E_LEAK_mV            = -54.387;
E_K_mV               = -88;
E_Na_mV              = 60;
G_LEAK               = 0.3;
G_K                  = 36;
G_Na                 = 120;

iInput_uA = @(t) 10 * (t > 100) - 10 * (t > 200) + 35 * (t > 300) - 35 * (t > 400);
simulationTime_ms    = 500;
simulationStepSize_ms = 0.0001;

commonArgs = {iInput_uA, REST_POTENTIAL_mV, E_LEAK_mV, E_K_mV, E_Na_mV, ...
              G_LEAK, G_K, G_Na, MEMBRANE_CAPACITANCE, ...
              simulationTime_ms, simulationStepSize_ms};

%% Run Both Solvers
tic;
% [vM_euler, n_euler, m_euler, h_euler, iLeak_euler, iK_euler, iNa_euler, t_euler] = ...
%     neuron_manual_integrate_faster(commonArgs{:});
[vM_euler, n_euler, m_euler, h_euler, iLeak_euler, iK_euler, iNa_euler, t_euler] = ...
    neuron_manual_integrate_faster(commonArgs{:});
time_euler = toc;

tic;
[vM_ode, n_ode, m_ode, h_ode, iLeak_ode, iK_ode, iNa_ode, t_ode] = ...
    neuron_variable_step_integrate(commonArgs{:});
time_ode = toc;

%% Compute Input Current & Conductances for Each Solver
inputCurrent_euler = iInput_uA(t_euler);
gK_euler   = G_K  * n_euler.^4;
gNa_euler  = G_Na * m_euler.^3 .* h_euler;

inputCurrent_ode = iInput_uA(t_ode);
gK_ode   = G_K  * n_ode.^4;
gNa_ode  = G_Na * m_ode.^3 .* h_ode;

%% Shared title suffix
eulerTitle = sprintf('Euler  |  T=%.0f ms, dt=%.4f ms  |  Runtime: %.3f s', ...
    simulationTime_ms, simulationStepSize_ms, time_euler);
odeTitle   = sprintf('ode15s  |  T=%.0f ms, dt=%.4f ms  |  Runtime: %.3f s', ...
    simulationTime_ms, simulationStepSize_ms, time_ode);

%% Helper
pad = @(y) [min(y(:)) - 0.05*(max(y(:))-min(y(:))+eps), ...
             max(y(:)) + 0.05*(max(y(:))-min(y(:))+eps)];

%% ---- Figure 1: Euler ------------------------------------------------
figure('Name', 'Euler Integration');
sgtitle(['Euler Integration  —  ' eulerTitle], 'Interpreter', 'none');

subplot(3,2,1);
plot(t_euler, vM_euler);
xlabel('Time (ms)'); ylabel('V_m (mV)');
title('Membrane Potential'); ylim('padded');

subplot(3,2,2);
plot(t_euler, inputCurrent_euler);
xlabel('Time (ms)'); ylabel('Current (\muA/cm^2)');
title('Injected Input Current'); ylim('padded');

subplot(3,2,3);
plot(t_euler, n_euler, 'b', t_euler, m_euler, 'r', t_euler, h_euler, 'g');
xlabel('Time (ms)'); ylabel('Gating value');
legend('n','m','h'); title('Gating Variables'); ylim('padded');

subplot(3,2,4);
plot(t_euler, iLeak_euler, '--', 'DisplayName', 'i_{Leak}'); hold on;
plot(t_euler, iK_euler,    '-.',  'DisplayName', 'i_{K}');
plot(t_euler, iNa_euler,   '-',   'DisplayName', 'i_{Na}'); hold off;
xlabel('Time (ms)'); ylabel('Current (\muA/cm^2)');
legend('Location','best'); title('Ionic Currents'); ylim('padded');

subplot(3,2,[5 6]);
plot(t_euler, gK_euler,  'b', 'DisplayName', 'g_K'); hold on;
plot(t_euler, gNa_euler, 'r', 'DisplayName', 'g_{Na}'); hold off;
xlabel('Time (ms)'); ylabel('Conductance (mS/cm^2)');
legend('Location','best'); title('Ionic Conductances'); ylim('padded');

%% ---- Figure 2: ode15s -----------------------------------------------
figure('Name', 'ode15s Integration');
sgtitle(['ode15s Integration  —  ' odeTitle], 'Interpreter', 'none');

subplot(3,2,1);
plot(t_ode, vM_ode);
xlabel('Time (ms)'); ylabel('V_m (mV)');
title('Membrane Potential'); ylim('padded');

subplot(3,2,2);
plot(t_ode, inputCurrent_ode);
xlabel('Time (ms)'); ylabel('Current (\muA/cm^2)');
title('Injected Input Current'); ylim('padded');

subplot(3,2,3);
plot(t_ode, n_ode, 'b', t_ode, m_ode, 'r', t_ode, h_ode, 'g');
xlabel('Time (ms)'); ylabel('Gating value');
legend('n','m','h'); title('Gating Variables'); ylim('padded');

subplot(3,2,4);
plot(t_ode, iLeak_ode, '--', 'DisplayName', 'i_{Leak}'); hold on;
plot(t_ode, iK_ode,    '-.',  'DisplayName', 'i_{K}');
plot(t_ode, iNa_ode,   '-',   'DisplayName', 'i_{Na}'); hold off;
xlabel('Time (ms)'); ylabel('Current (\muA/cm^2)');
legend('Location','best'); title('Ionic Currents'); ylim('padded');

subplot(3,2,[5 6]);
plot(t_ode, gK_ode,  'b', 'DisplayName', 'g_K'); hold on;
plot(t_ode, gNa_ode, 'r', 'DisplayName', 'g_{Na}'); hold off;
xlabel('Time (ms)'); ylabel('Conductance (mS/cm^2)');
legend('Location','best'); title('Ionic Conductances'); ylim('padded');

%% ---- Figure 3: Error between solvers --------------------------------
% Interpolate ode15s onto Euler time grid for a point-wise comparison
vM_ode_interp    = interp1(t_ode, vM_ode,    t_euler, 'linear');
n_ode_interp     = interp1(t_ode, n_ode,     t_euler, 'linear');
m_ode_interp     = interp1(t_ode, m_ode,     t_euler, 'linear');
h_ode_interp     = interp1(t_ode, h_ode,     t_euler, 'linear');
iK_ode_interp    = interp1(t_ode, iK_ode,    t_euler, 'linear');
iNa_ode_interp   = interp1(t_ode, iNa_ode,   t_euler, 'linear');
gK_ode_interp    = interp1(t_ode, gK_ode,    t_euler, 'linear');
gNa_ode_interp   = interp1(t_ode, gNa_ode,   t_euler, 'linear');

figure('Name', 'Solver Error');
sgtitle('Absolute Error: Euler vs ode15s', 'Interpreter', 'none');

subplot(3,2,1);
plot(t_euler, abs(vM_euler - vM_ode_interp));
xlabel('Time (ms)'); ylabel('\DeltaV_m (mV)');
title('Membrane Potential Error'); ylim('padded');

subplot(3,2,2);
plot(t_euler, abs(n_euler - n_ode_interp), 'b', ...
     t_euler, abs(m_euler - m_ode_interp), 'r', ...
     t_euler, abs(h_euler - h_ode_interp), 'g');
xlabel('Time (ms)'); ylabel('\Delta Gating');
legend('n','m','h'); title('Gating Variable Error'); ylim('padded');

subplot(3,2,3);
plot(t_euler, abs(iK_euler - iK_ode_interp), 'b', ...
     t_euler, abs(iNa_euler - iNa_ode_interp), 'r');
xlabel('Time (ms)'); ylabel('\Delta Current (\muA/cm^2)');
legend('i_K','i_{Na}'); title('Ionic Current Error'); ylim('padded');

subplot(3,2,4);
plot(t_euler, abs(gK_euler - gK_ode_interp), 'b', ...
     t_euler, abs(gNa_euler - gNa_ode_interp), 'r');
xlabel('Time (ms)'); ylabel('\Delta Conductance (mS/cm^2)');
legend('g_K','g_{Na}'); title('Conductance Error'); ylim('padded');

subplot(3,2,[5 6]);
rmsVm = sqrt(mean((vM_euler - vM_ode_interp).^2));
rmsN  = sqrt(mean((n_euler  - n_ode_interp ).^2));
rmsM  = sqrt(mean((m_euler  - m_ode_interp ).^2));
rmsH  = sqrt(mean((h_euler  - h_ode_interp ).^2));
bar(categorical({'V_m','n','m','h'}), [rmsVm, rmsN, rmsM, rmsH]);
ylabel('RMS Error'); title('RMS Error Summary');