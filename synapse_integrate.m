%synapse_integration
pulseAmplitude = 30;     % mV
pulseStart = 0;          % ms
pulseWidth = 5;          % ms

V_pre = @(t) pulseAmplitude .* (t >= pulseStart & t < (pulseStart + pulseWidth));

% Input spike train driving h(t)
h_in = @(t) V_pre(t);


% Initialize parameters for synapse integration
timeStepSize = 0.01; % Time step for simulation
simulationTime = 75;   % Total simulation time

g_bar_syn_s = 0.00006; %mS
tau_rise = 2;  %ms 
tau_decay = 11; %ms (weird units but avoid an additional mult)
E_syn = -62.5; % (mV)

t_peak = tau_decay * tau_rise * log(tau_decay / tau_rise)/(tau_decay - tau_rise);

f_norm = 1/((-exp(-t_peak / tau_rise) + exp(-t_peak/tau_decay)));


%This part is necessary for spike-mediated synapses between oscillator
%neurons

% v should be v_pre

function dfdt = odesys(t, state, tau_rise, tau_decay, h_in)
    
    h_ = state(1);
    g_ = state(2);

    dh = -h_ / tau_rise + h_in(t);
    dg = -g_ / tau_decay + h_;

    dfdt = [dh; dg];

end

h_0 = 0;
g_0 = 0;

state0 = [h_0,g_0];

tSpan = 0: timeStepSize: simulationTime;
options = odeset('RelTol',1e-6, 'AbsTol', 1e-8);
% Set up the simulation using ode45
% Pass parameters into odesys using anonymous function
[timeStamps, solution] = ode15s(@(t,state) odesys(t, state, tau_rise, tau_decay, h_in), [0 simulationTime], state0,    options);

timeStamps = timeStamps';
h = solution(:, 1)';
g = solution(:, 2)';

%Calculate results of g
g_syn = g_bar_syn_s * f_norm * g;

%uA
Isyn = g_syn.*(V_pre(timeStamps) - E_syn);

figure;
plot(timeStamps, Isyn)
xlabel('Time (ms)')
ylabel('I_{syn} (\muA)')
title('Synaptic Current')
grid on

figure;
plot(timeStamps, g_syn)
xlabel('Time (ms)')
ylabel('g_{syn} (mS)')
title('Synaptic Conductance')
grid on