clc;
clear;
%% Neuron Integration Runner
MEMBRANE_CAPACITANCE = 1; % uF/cm^2
REST_POTENTIAL_mV = -65;
E_LEAK_mV = -54.387;
E_K_mV = -88;
E_Na_mV = 60;
G_LEAK = 0.3;
G_K = 36;
G_Na = 120;

iInput_uA = @(t) 10 * (t > 100) - 10 * (t > 200) + 35 * (t > 300) - 35 * (t > 400);
%iInput_uA = @(t) (15*cos(2*pi*(1/50)*t) + 15) .* (t > 50);

simulationTime_ms = 500;
simulationStepSize_ms = 0.005;


%[vMembrane,n,m,h,iLeak,iK,iNa,timeStamps] = neuron_manual_integrate(iInput_uA,REST_POTENTIAL_mV,E_LEAK_mV,E_K_mV,E_Na_mV,G_LEAK,G_K,G_Na,MEMBRANE_CAPACITANCE,simulationTime_ms,simulationStepSize_ms);
[vMembrane,n,m,h,iLeak,iK,iNa,timeStamps] = neuron_variable_step_integrate(iInput_uA,REST_POTENTIAL_mV,E_LEAK_mV,E_K_mV,E_Na_mV,G_LEAK,G_K,G_Na,MEMBRANE_CAPACITANCE,simulationTime_ms,simulationStepSize_ms);



% Compute input current trace for plotting
inputCurrent = iInput_uA(timeStamps);

% Create figure arranged as 2x2
figure;

% Helper to add small margins to y-limits
add_margin = @(y) [min(y(:)) - 0.05*(max(y(:))-min(y(:)) + eps), ...
                   max(y(:)) + 0.05*(max(y(:))-min(y(:)) + eps)];
% Plot membrane potential (top-left)
subplot(2,2,1);
plot(timeStamps, vMembrane);
xlabel('Time (ms)');
ylabel('V_m (mV)');
title('Membrane Potential');
ylim('padded');


% Plot input current (top-right)
subplot(2,2,2);
plot(timeStamps, inputCurrent);
xlabel('Time (ms)');
ylabel('Input Current (\muA/cm^2)');
title('Injected Input Current');
ylim('padded');


% Plot gating variables n, m, h (bottom-left)
subplot(2,2,3);
plot(timeStamps, n, 'b', timeStamps, m, 'r', timeStamps, h, 'g');
xlabel('Time (ms)');
ylabel('Gating value');
legend('n','m','h');
title('Gating Variables');
ylim('padded');


% Plot currents iLeak, iK, iNa together (bottom-right)
subplot(2,2,4);
plot(timeStamps, iLeak, '--','DisplayName','i_{Leak}');
hold on;
plot(timeStamps, iK, '-.','DisplayName','i_{K}');
plot(timeStamps, iNa, '-','DisplayName','i_{Na}');
hold off;
xlabel('Time (ms)');
ylabel('Current (\muA/cm^2)');
legend('Location','best');
title('Ionic Currents');
ylim('padded');
