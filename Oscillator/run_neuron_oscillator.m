%All voltages in V, all conductances in S, all times in seconds, capacitances in F 

restPotential = -0.060;

eqPotLeak = restPotential;  gLeak     = 8e-9;
eqPotNa   =  0.045;  gNaFast   = 200e-9;  gNaPers   = 7e-9;
eqPotK    = -0.070;  gKDelRect = 100e-9;  gKPers    = 80e-9;  gKFast = 80e-9;
eqPotCa   =  0.135;  gCaFast   = 5e-9;    gCaSlow   = 3.2e-9;
eqPotHypPol = -0.021; gHypPol  = 4e-9;
eqPotSyn  = -0.0625; gSynS     = 60e-9;

tau_rise  = 0.002;   tau_decay = 0.011; %(seconds) 
h0Thresh  = -0.020;
capacitance   = 5e-10;
simulationTime = 30;

% ---- Corrected call (all SI: V, S, F, seconds) ----
tic

[timeStamps, solution] = neuron_oscillator_variable_step_integrate( ...
    restPotential, ...
    eqPotLeak, gLeak, ...
    eqPotNa, gNaFast, gNaPers, ...
    eqPotK, gKDelRect, gKPers, gKFast, ...
    eqPotCa, gCaFast, gCaSlow, ...
    eqPotHypPol, gHypPol, ...
    eqPotSyn, gSynS, ...
    tau_rise, tau_decay, h0Thresh, ...
    capacitance, simulationTime);
runtime = toc;

%% Output
timeStamps  = timeStamps';
L_vMembrane = solution(:, 1)';
R_vMembrane = solution(:, 2)';
L_mNaFast   = solution(:, 3)';
R_mNaFast   = solution(:, 4)';
L_hNaFast   = solution(:, 5)';
R_hNaFast   = solution(:, 6)';
L_mNaPers   = solution(:, 7)';
R_mNaPers   = solution(:, 8)';
L_mCaFast   = solution(:, 9)';
R_mCaFast   = solution(:,10)';
L_hCaFast   = solution(:,11)';
R_hCaFast   = solution(:,12)';
L_mCaSlow   = solution(:,13)';
R_mCaSlow   = solution(:,14)';
L_hCaSlow   = solution(:,15)';
R_hCaSlow   = solution(:,16)';
L_mHypPol   = solution(:,17)';
R_mHypPol   = solution(:,18)';
L_mKDelRect = solution(:,19)';
R_mKDelRect = solution(:,20)';
L_hKDelRect = solution(:,21)';
R_hKDelRect = solution(:,22)';
L_mKPers    = solution(:,23)';
R_mKPers    = solution(:,24)';
L_mKFast    = solution(:,25)';
R_mKFast    = solution(:,26)';
L_hKFast    = solution(:,27)';
R_hKFast    = solution(:,28)';
L_syn_h     = solution(:,29)';
R_syn_h     = solution(:,30)';
L_syn_g     = solution(:,31)';
R_syn_g     = solution(:,32)';

% Compute currents time-series for each timestep using previously extracted state variables
nT = numel(timeStamps);

% Preallocate current arrays
L_iNaFast  = zeros(1,nT); R_iNaFast  = zeros(1,nT);
L_iNaPers  = zeros(1,nT); R_iNaPers  = zeros(1,nT);
L_iCaFast  = zeros(1,nT); R_iCaFast  = zeros(1,nT);
L_iCaSlow  = zeros(1,nT); R_iCaSlow  = zeros(1,nT);
L_iHypPol  = zeros(1,nT); R_iHypPol  = zeros(1,nT);
L_iKDelRect = zeros(1,nT); R_iKDelRect = zeros(1,nT);
L_iKPers   = zeros(1,nT); R_iKPers   = zeros(1,nT);
L_iKFast   = zeros(1,nT); R_iKFast   = zeros(1,nT);
L_iLeak    = zeros(1,nT); R_iLeak    = zeros(1,nT);
L_iSynS    = zeros(1,nT); R_iSynS    = zeros(1,nT);


% ---- Recompute syn_f_norm (mirrors logic inside integrator) ----
syn_t_peak = tau_decay * tau_rise * log(tau_decay / tau_rise) / (tau_decay - tau_rise);
syn_f_norm = 1 / (-exp(-syn_t_peak / tau_rise) + exp(-syn_t_peak / tau_decay));


% Loop vectorized across timesteps (vector operations)
% Use elementwise operations on vectors already available
% all of them are in Amps
L_iNaFast  = gNaFast .* (L_mNaFast .^ 3) .* L_hNaFast .* (L_vMembrane - eqPotNa);
R_iNaFast  = gNaFast .* (R_mNaFast .^ 3) .* R_hNaFast .* (R_vMembrane - eqPotNa);

L_iNaPers  = gNaPers .* L_mNaPers .* (L_vMembrane - eqPotNa);
R_iNaPers  = gNaPers .* R_mNaPers .* (R_vMembrane - eqPotNa);

L_iCaFast  = gCaFast .* (L_mCaFast .^ 2) .* L_hCaFast .* (L_vMembrane - eqPotCa);
R_iCaFast  = gCaFast .* (R_mCaFast .^ 2) .* R_hCaFast .* (R_vMembrane - eqPotCa);

L_iCaSlow  = gCaSlow .* (L_mCaSlow .^ 2) .* L_hCaSlow .* (L_vMembrane - eqPotCa);
R_iCaSlow  = gCaSlow .* (R_mCaSlow .^ 2) .* R_hCaSlow .* (R_vMembrane - eqPotCa);

L_iHypPol  = gHypPol .* (L_mHypPol .^ 2) .* (L_vMembrane - eqPotHypPol);
R_iHypPol  = gHypPol .* (R_mHypPol .^ 2) .* (R_vMembrane - eqPotHypPol);

L_iKDelRect = gKDelRect .* (L_mKDelRect .^ 2) .* L_hKDelRect .* (L_vMembrane - eqPotK);
R_iKDelRect = gKDelRect .* (R_mKDelRect .^ 2) .* R_hKDelRect .* (R_vMembrane - eqPotK);

L_iKPers   = gKPers .* (L_mKPers .^ 2) .* (L_vMembrane - eqPotK);
R_iKPers   = gKPers .* (R_mKPers .^ 2) .* (R_vMembrane - eqPotK);

L_iKFast   = gKFast .* (L_mKFast .^ 2) .* L_hKFast .* (L_vMembrane - eqPotK);
R_iKFast   = gKFast .* (R_mKFast .^ 2) .* R_hKFast .* (R_vMembrane - eqPotK);

L_iLeak    = gLeak .* (L_vMembrane - eqPotLeak);
R_iLeak    = gLeak .* (R_vMembrane - eqPotLeak);

L_iSynS    = gSynS .* syn_f_norm .* L_syn_g .* (L_vMembrane - eqPotSyn);
R_iSynS    = gSynS .* syn_f_norm .* R_syn_g .* (R_vMembrane - eqPotSyn);

% Optionally collect into a struct for convenience
currents.time = timeStamps;
currents.L_iNaFast  = L_iNaFast;  currents.R_iNaFast  = R_iNaFast;
currents.L_iNaPers  = L_iNaPers;  currents.R_iNaPers  = R_iNaPers;
currents.L_iCaFast  = L_iCaFast;  currents.R_iCaFast  = R_iCaFast;
currents.L_iCaSlow  = L_iCaSlow;  currents.R_iCaSlow  = R_iCaSlow;
currents.L_iHypPol  = L_iHypPol;  currents.R_iHypPol  = R_iHypPol;
currents.L_iKDelRect = L_iKDelRect; currents.R_iKDelRect = R_iKDelRect;
currents.L_iKPers   = L_iKPers;   currents.R_iKPers   = R_iKPers;
currents.L_iKFast   = L_iKFast;   currents.R_iKFast   = R_iKFast;
currents.L_iLeak    = L_iLeak;    currents.R_iLeak    = R_iLeak;
currents.L_iSynS    = L_iSynS;    currents.R_iSynS    = R_iSynS;



% ---- Plot (output is in V, convert to mV) ----
figure;
title('Membrane Potentials')
subplot(2,1,1);
plot(currents.time, L_vMembrane * 1e3, 'b', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
title('Left Neuron'); grid on;

subplot(2,1,2);
plot(currents.time, R_vMembrane * 1e3, 'r', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Membrane Potential (mV)');
title('Right Neuron'); grid on;

% ---- Additional plots: conductances and currents grouped by figure ----
to_nS = 1e9;
to_nA = 1e9;
t = currents.time;

% Include runtime in the first figure title
runtimeStr = sprintf(' (runtime: %.3f s)', runtime);

% Helper to create a 2-row grid of subplots for a set of traces on a given figure
function create_pair_subplot(figHandle, idx, t, gL, gR, gLabel, iL, iR, iLabel, titleSuffix)
    figure(figHandle);
    ax1 = subplot(2,1,1);
    plot(t, gL, 'b', 'LineWidth',1.2); hold on;
    plot(t, gR, 'r', 'LineWidth',1.2); hold off;
    xlabel('Time (s)'); ylabel(gLabel); title([titleSuffix ' - Conductance']); grid on;
    legend('L','R');
    ax2 = subplot(2,1,2);
    plot(t, iL, 'b', 'LineWidth',1.2); hold on;
    plot(t, iR, 'r', 'LineWidth',1.2); hold off;
    xlabel('Time (s)'); ylabel(iLabel); title([titleSuffix ' - Current']); grid on;
    legend('L','R');
end

% --- Sodium: put fast and persistent in one figure with two columns (each column is a pair of subplots) ---
figNa = figure('Name', ['Sodium (Na)' runtimeStr],'NumberTitle','off');
% Left column: Fast sodium
subplot(2,2,1);
plot(t, gNaFast .* (L_mNaFast .^ 3) .* L_hNaFast * to_nS, 'b', 'LineWidth',1.2); hold on;
plot(t, gNaFast .* (R_mNaFast .^ 3) .* R_hNaFast * to_nS, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('g_{Na,fast} (nS)'); title('Na - Fast Conductance'); grid on; legend('L','R');

subplot(2,2,3);
plot(t, L_iNaFast*to_nA, 'b', 'LineWidth',1.2); hold on;
plot(t, R_iNaFast*to_nA, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('I_{Na,fast} (nA)'); title('Na - Fast Current'); grid on; legend('L','R');

% Right column: Persistent sodium
subplot(2,2,2);
plot(t, gNaPers .* L_mNaPers * to_nS, 'b', 'LineWidth',1.2); hold on;
plot(t, gNaPers .* R_mNaPers * to_nS, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('g_{Na,pers} (nS)'); title('Na - Persistent Conductance'); grid on; legend('L','R');

subplot(2,2,4);
plot(t, L_iNaPers*to_nA, 'b', 'LineWidth',1.2); hold on;
plot(t, R_iNaPers*to_nA, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('I_{Na,pers} (nA)'); title('Na - Persistent Current'); grid on; legend('L','R');

% --- Calcium: fast and slow in one figure ---
figCa = figure('Name','Calcium (Ca)','NumberTitle','off');

% Fast Ca
subplot(2,2,1);
plot(t, gCaFast .* (L_mCaFast .^ 2) .* L_hCaFast * to_nS, 'b', 'LineWidth',1.2); hold on;
plot(t, gCaFast .* (R_mCaFast .^ 2) .* R_hCaFast * to_nS, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('g_{Ca,fast} (nS)'); title('Ca - Fast Conductance'); grid on; legend('L','R');

subplot(2,2,3);
plot(t, L_iCaFast*to_nA, 'b', 'LineWidth',1.2); hold on;
plot(t, R_iCaFast*to_nA, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('I_{Ca,fast} (nA)'); title('Ca - Fast Current'); grid on; legend('L','R');

% Slow Ca
subplot(2,2,2);
plot(t, gCaSlow .* (L_mCaSlow .^ 2) .* L_hCaSlow * to_nS, 'b', 'LineWidth',1.2); hold on;
plot(t, gCaSlow .* (R_mCaSlow .^ 2) .* R_hCaSlow * to_nS, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('g_{Ca,slow} (nS)'); title('Ca - Slow Conductance'); grid on; legend('L','R');

subplot(2,2,4);
plot(t, L_iCaSlow*to_nA, 'b', 'LineWidth',1.2); hold on;
plot(t, R_iCaSlow*to_nA, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('I_{Ca,slow} (nA)'); title('Ca - Slow Current'); grid on; legend('L','R');

% --- Potassium: delayed-rectifier, persistent, fast in one figure (3 columns) ---
figK = figure('Name','Potassium (K)','NumberTitle','off');

% Delayed-rectifier (column 1)
subplot(2,3,1);
plot(t, gKDelRect .* (L_mKDelRect .^ 2) .* L_hKDelRect * to_nS, 'b', 'LineWidth',1.2); hold on;
plot(t, gKDelRect .* (R_mKDelRect .^ 2) .* R_hKDelRect * to_nS, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('g_{K,del} (nS)'); title('K - DelayedRect Conductance'); grid on; legend('L','R');

subplot(2,3,4);
plot(t, L_iKDelRect*to_nA, 'b', 'LineWidth',1.2); hold on;
plot(t, R_iKDelRect*to_nA, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('I_{K,del} (nA)'); title('K - DelayedRect Current'); grid on; legend('L','R');

% Persistent (column 2)
subplot(2,3,2);
plot(t, gKPers .* (L_mKPers .^ 2) * to_nS, 'b', 'LineWidth',1.2); hold on;
plot(t, gKPers .* (R_mKPers .^ 2) * to_nS, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('g_{K,pers} (nS)'); title('K - Persistent Conductance'); grid on; legend('L','R');

subplot(2,3,5);
plot(t, L_iKPers*to_nA, 'b', 'LineWidth',1.2); hold on;
plot(t, R_iKPers*to_nA, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('I_{K,pers} (nA)'); title('K - Persistent Current'); grid on; legend('L','R');

% Fast (column 3)
subplot(2,3,3);
plot(t, gKFast .* (L_mKFast .^ 2) .* L_hKFast * to_nS, 'b', 'LineWidth',1.2); hold on;
plot(t, gKFast .* (R_mKFast .^ 2) .* R_hKFast * to_nS, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('g_{K,fast} (nS)'); title('K - Fast Conductance'); grid on; legend('L','R');

subplot(2,3,6);
plot(t, L_iKFast*to_nA, 'b', 'LineWidth',1.2); hold on;
plot(t, R_iKFast*to_nA, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('I_{K,fast} (nA)'); title('K - Fast Current'); grid on; legend('L','R');

% --- Leak and Hyperpolarizing on the same figure (two rows) ---
figMisc = figure('Name','Leak & Hyperpolarizing','NumberTitle','off');

% Leak (row 1)
subplot(2,1,1);
plot(t, repmat(gLeak*to_nS,1,nT), 'k--', 'LineWidth',1.2); hold on;
plot(t, repmat(gLeak*to_nS,1,nT), 'k--', 'LineWidth',1.2); hold off; % identical for L and R
xlabel('Time (s)'); ylabel('g_{leak} (nS)'); title('Leak Conductance'); grid on; legend('L','R');
subplot(2,1,2);
plot(t, L_iLeak*to_nA, 'b', 'LineWidth',1.2); hold on;
plot(t, R_iLeak*to_nA, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('I_{leak} (nA)'); title('Leak Current'); grid on; legend('L','R');

% New figure for Hyperpolarizing plotted beneath leak's current for clarity (reuse same figMisc as subplot positions)
figHyper = figure('Name','Hyperpolarizing','NumberTitle','off');
subplot(2,1,1);
plot(t, (gHypPol .* (L_mHypPol .^ 2)) * to_nS, 'b', 'LineWidth',1.2); hold on;
plot(t, (gHypPol .* (R_mHypPol .^ 2)) * to_nS, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('g_{hyp} (nS)'); title('Hyperpolarizing Conductance'); grid on; legend('L','R');

subplot(2,1,2);
plot(t, L_iHypPol*to_nA, 'b', 'LineWidth',1.2); hold on;
plot(t, R_iHypPol*to_nA, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('I_{hyp} (nA)'); title('Hyperpolarizing Current'); grid on; legend('L','R');

% --- Synapse in its own figure ---
figSyn = figure('Name','Synapse','NumberTitle','off');
subplot(2,1,1);
plot(t, L_syn_g * to_nS, 'b', 'LineWidth',1.2); hold on;
plot(t, R_syn_g * to_nS, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('g_{syn} (nS)'); title('Synapse Conductance'); grid on; legend('L','R');

subplot(2,1,2);
plot(t, L_iSynS*to_nA, 'b', 'LineWidth',1.2); hold on;
plot(t, R_iSynS*to_nA, 'r', 'LineWidth',1.2); hold off;
xlabel('Time (s)'); ylabel('I_{syn} (nA)'); title('Synapse Current'); grid on; legend('L','R');