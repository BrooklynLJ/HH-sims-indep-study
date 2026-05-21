function [vMembrane, n, m, h, iLeak, iK, iNa, timeStamps] = neuron_oscillator_variable_step_integrate(...
    iInput, restPotential, ...
    eqPotLeak, gLeak, ...
    eqPotNa, gNaFast, gNaPers, ...
    eqPotK, gKDelRect, gKPers, gKFast, ...
    eqPotCa, gCaFast, gCaSlow, ...
    eqPotHypPol,  gHypPol, ...
    eqPotSyn, gSynS, ...
    tau_decay, tau_rise, ...
    capacitance, ...
    simulationTime, timeStepSize)

% Description: ODE15s integration of the Hodgkin-Huxley model.
% Units: Time (ms), Current (uA/cm^2), Potential (mV), Conductance (mS/cm^2)
% Note: timeStepSize defines the output evaluation points, not the solver step size.

%% 1. Simulation Setup
% totalNumSteps = simulationTime / timeStepSize;
% dt = timeStepSize;

%% 2. Gating Variable Dynamics (Hodgkin-Huxley)
fInf = @(a,b,V) 1/(1+exp(a*(V+b)));
fhInf = @(V) 1/(1 + 2*exp(180*(V+0.049)) + exp(500*(V+0.049)));
tau = @(a, b, c, d, V) c + d /(1 + exp(a*(V+b)));
%tauhKF = @(V) 1.5 + 8/1+exp(-100 *(V+0.022)) + (-2.2)/cosh(100 *(V+0.04));

%Set up all DE auxillary functions
%Na Fast
mNaFast_inf = @(V) fInf(-150, 0.029, V);
taumNaFast = @(V) 0.0001;

hNaFast_inf = @(V) fInf(500, 0.03, V);
tauhNaFast = @(V) 0.004 + 0.006/(1 + exp(500 *(V + 0.028))) + 0.01/(cosh(300 * (V + 0.027)));

%Na Persistent 
mNaPers_inf = @(V) fInf(-120, 0.039, V);
taumNaPers = @(V) tau(400, 0.057, 0.01, 0.2, V);

% Ca Fast
mCaFast_inf = @(V) fInf(-600, 0.0467, V);
taumCaFast = @(V) 0.011 + 0.024/(cosh(-330*(V+0.0467)));

hCaFast_inf = @(V) fInf(350, 0.0555, V);
tauhCaFast = @(V) tau(270, 0.055, 0.06, 0.31, V);

% Ca slow
mCaSlow_inf = @(V) fInf(-420, 0.0472, V);
taumCaSlow = @(V) tau(-400,0.0487, 0.005, 0.134, V);

hCaSlow_inf = @(V) fInf(360, 0.055, V);
tauhCaSlow = @(V) tau(-250, 0.043, 0.2, 5.25, V);

% Hyperpolarization activated cation 
mHypPol_inf = @(V) fhInf(V);
taumHypPol = @(V) tau(-100, 0.073, 0.7, 1.7, V);

% Delayed Rectifier Like K
mKDelRect_inf = @(V) fInf(-143, 0.021, V);
taumKDelRect = @(V) tau(150, 0.016, 0.001, 0.011, V);

hKDelRect_inf = @(V) fInf(111,0.028,V);
tauhKDelRect = @(V) tau(-143,0.013, 0.5, 0.2, V);

% persistent Potassium Current K2;
mKPers_inf = @(V) fInf(-83,0.02, V);
taumKPers = @(V) tau(200, 0.035, 0.057, 0.043, V);

% Fast potasium current
mKFast_inf = @(V) fInf(-130,0.044,V);
taumKFast_inf = @(V) tau(200, 0.03, 0.005, 0.011, V);

hKFast_inf = @(V) fInf(160, 0.063, V);
tauhKFast = @(V) tau(-300, 0.055, 0.026, 0.0085, V);


%% Synapse model

syn_t_peak = tau_decay * tau_rise * log(tau_decay / tau_rise)/(tau_decay - tau_rise);
syn_f_norm = 1/((-exp(-t_peak / tau_rise) + exp(-t_peak/tau_decay)));




%% Freaky matlab ode function declaration:
    function dfdt = odesys(t,state)
        L_vMembrane = state(1);
        R_vMembrane = state(2);
        
        L_mNaFast = state(3);
        d_L_mNaFast = (mNaFast_inf(L_vMembrane) - L_mNaFast) / taumNaFast(L_vMembrane);
        
        R_mNaFast = state(4);
        d_R_mNaFast = (mNaFast_inf(R_vMembrane) - R_mNaFast) / taumNaFast(R_vMembrane);
        
        L_hNaFast = state(5);
        d_L_hNaFast = (hNaFast_inf(L_vMembrane) - L_hNaFast) / tauhNaFast(L_vMembrane);
        
        R_hNaFast = state(6);
        d_R_hNaFast = (hNaFast_inf(R_vMembrane) - R_hNaFast) / tauhNaFast(R_vMembrane);
        
        % Na Persistent
        L_mNaPers = state(7);
        d_L_mNaPers = (mNaPers_inf(L_vMembrane) - L_mNaPers) / taumNaPers(L_vMembrane);
        
        R_mNaPers = state(8);
        d_R_mNaPers = (mNaPers_inf(R_vMembrane) - R_mNaPers) / taumNaPers(R_vMembrane);
        
        % Ca Fast
        L_mCaFast = state(9);
        d_L_mCaFast = (mCaFast_inf(L_vMembrane) - L_mCaFast) / taumCaFast(L_vMembrane);
        
        R_mCaFast = state(10);
        d_R_mCaFast = (mCaFast_inf(R_vMembrane) - R_mCaFast) / taumCaFast(R_vMembrane);
        
        L_hCaFast = state(11);
        d_L_hCaFast = (hCaFast_inf(L_vMembrane) - L_hCaFast) / tauhCaFast(L_vMembrane);
        
        R_hCaFast = state(12);
        d_R_hCaFast = (hCaFast_inf(R_vMembrane) - R_hCaFast) / tauhCaFast(R_vMembrane);
        
        % Ca slow
        L_mCaSlow = state(13);
        d_L_mCaSlow = (mCaSlow_inf(L_vMembrane) - L_mCaSlow) / taumCaSlow(L_vMembrane);
        
        R_mCaSlow = state(14);
        d_R_mCaSlow = (mCaSlow_inf(R_vMembrane) - R_mCaSlow) / taumCaSlow(R_vMembrane);
        
        L_hCaSlow = state(15);
        d_L_hCaSlow = (hCaSlow_inf(L_vMembrane) - L_hCaSlow) / tauhCaSlow(L_vMembrane);
        
        R_hCaSlow = state(16);
        d_R_hCaSlow = (hCaSlow_inf(R_vMembrane) - R_hCaSlow) / tauhCaSlow(R_vMembrane);
        
        % Hyperpolarization activated cation
        L_mHypPol = state(17);
        d_L_mHypPol = (mHypPol_inf(L_vMembrane) - L_mHypPol) / taumHypPol(L_vMembrane);
        
        R_mHypPol = state(18);
        d_R_mHypPol = (mHypPol_inf(R_vMembrane) - R_mHypPol) / taumHypPol(R_vMembrane);
        
        % Delayed Rectifier Like K
        L_mKDelRect = state(19);
        d_L_mKDelRect = (mKDelRect_inf(L_vMembrane) - L_mKDelRect) / taumKDelRect(L_vMembrane);
        
        R_mKDelRect = state(20);
        d_R_mKDelRect = (mKDelRect_inf(R_vMembrane) - R_mKDelRect) / taumKDelRect(R_vMembrane);
        
        L_hKDelRect = state(21);
        d_L_hKDelRect = (hKDelRect_inf(L_vMembrane) - L_hKDelRect) / tauhKDelRect(L_vMembrane);
        
        R_hKDelRect = state(22);
        d_R_hKDelRect = (hKDelRect_inf(R_vMembrane) - R_hKDelRect) / tauhKDelRect(R_vMembrane);
        
        % persistent Potassium Current K2
        L_mKPers = state(23);
        d_L_mKPers = (mKPers_inf(L_vMembrane) - L_mKPers) / taumKPers(L_vMembrane);
        
        R_mKPers = state(24);
        d_R_mKPers = (mKPers_inf(R_vMembrane) - R_mKPers) / taumKPers(R_vMembrane);
        
        % Fast potassium current
        L_mKFast = state(25);
        d_L_mKFast = (mKFast_inf(L_vMembrane) - L_mKFast) / taumKFast_inf(L_vMembrane);
        
        R_mKFast = state(26);
        d_R_mKFast = (mKFast_inf(R_vMembrane) - R_mKFast) / taumKFast_inf(R_vMembrane);
        
        L_hKFast = state(27);
        d_L_hKFast = (hKFast_inf(L_vMembrane) - L_hKFast) / tauhKFast(L_vMembrane);
        
        R_hKFast = state(28);
        d_R_hKFast = (hKFast_inf(R_vMembrane) - R_hKFast) / tauhKFast(R_vMembrane);

        % Synapse stuff
        L_syn_h = state(29);
        d_L_syn_h = - L_syn_h / tau_rise + R_vMembrane;

        R_syn_h = state(30);
        d_R_syn_h = - R_syn_h / tau_rise + L_vMembrane;
        
        L_syn_g = state(31);
        d_L_syn_g = -L_syn_g / tau_decay + L_syn_h;

        R_syn_g = state(32);
        d_R_syn_g = - R_syn_g / tau_decay + R_syn_h;


% --- Left and Right Currents ---
        
        % Fast Sodium
        L_iNaFast = gNaFast * (L_mNaFast ^ 3) * L_hNaFast * (L_vMembrane - eqPotNa);
        R_iNaFast = gNaFast * (R_mNaFast ^ 3) * R_hNaFast * (R_vMembrane - eqPotNa);
        
        % Persistent Sodium
        L_iNaPers = gNaPers * L_mNaPers * (L_vMembrane - eqPotNa);
        R_iNaPers = gNaPers * R_mNaPers * (R_vMembrane - eqPotNa);
        
        % Fast Calcium
        L_iCaFast = gCaFast * (L_mCaFast ^ 2) * L_hCaFast * (L_vMembrane - eqPotCa);
        R_iCaFast = gCaFast * (R_mCaFast ^ 2) * R_hCaFast * (R_vMembrane - eqPotCa);
        
        % Slow Calcium
        L_iCaSlow = gCaSlow * (L_mCaSlow ^ 2) * L_hCaSlow * (L_vMembrane - eqPotCa);
        R_iCaSlow = gCaSlow * (R_mCaSlow ^ 2) * R_hCaSlow * (R_vMembrane - eqPotCa);
        
        % Hyperpolarization Activated Cation
        L_iHypPol = gHypPol * (L_mHypPol ^ 2) * (L_vMembrane - eqPotHypPol);
        R_iHypPol = gHypPol * (R_mHypPol ^ 2) * (R_vMembrane - eqPotHypPol);
        
        % Delayed Rectifier Like K
        L_iKDelRect = gKDelRect * (L_mKDelRect ^ 2) * L_hKDelRect * (L_vMembrane - eqPotK);
        R_iKDelRect = gKDelRect * (R_mKDelRect ^ 2) * R_hKDelRect * (R_vMembrane - eqPotK);
        
        % Persistent Potassium Current K2
        L_iKPers = gKPers * (L_mKPers ^ 2) * (L_vMembrane - eqPotK);
        R_iKPers = gKPers * (R_mKPers ^ 2) * (R_vMembrane - eqPotK);
        
        % Fast Potassium
        L_iKFast = gKFast * (L_mKFast ^ 2) * L_hKFast * (L_vMembrane - eqPotK);
        R_iKFast = gKFast * (R_mKFast ^ 2) * R_hKFast * (R_vMembrane - eqPotK);
        
        % Leak
        L_iLeak = gLeak * (L_vMembrane - eqPotLeak);
        R_iLeak = gLeak * (R_vMembrane - eqPotLeak);

        L_iSynS = gSynS * syn_f_norm * L_syn_g * (R_vMembrane - eqPotSyn);
        R_iSynS = gSynS * syn_f_norm * R_syn_g * (L_vMembrane - eqPotSyn);
        

        d_L_vMembrane = - (1 / capacitance) * (L_iNaFast + L_iNaPers + L_iCaFast + L_iCaSlow + L_iHypPol + L_iKDelRect + L_iKPers + L_iKFast + L_iLeak + L_iSynS);
        d_R_vMembrane = - (1 / capacitance) * (R_iNaFast + R_iNaPers + R_iCaFast + R_iCaSlow + R_iHypPol + R_iKDelRect + R_iKPers + R_iKFast + R_iLeak + R_iSynS);

dfdt = [d_L_vMembrane; d_R_vMembrane; ...
                d_L_mNaFast;   d_R_mNaFast;   ...
                d_L_hNaFast;   d_R_hNaFast;   ...
                d_L_mNaPers;   d_R_mNaPers;   ...
                d_L_mCaFast;   d_R_mCaFast;   ...
                d_L_hCaFast;   d_R_hCaFast;   ...
                d_L_mCaSlow;   d_R_mCaSlow;   ...
                d_L_hCaSlow;   d_R_hCaSlow;   ...
                d_L_mHypPol;   d_R_mHypPol;   ...
                d_L_mKDelRect; d_R_mKDelRect; ...
                d_L_hKDelRect; d_R_hKDelRect; ...
                d_L_mKPers;    d_R_mKPers;    ...
                d_L_mKFast;    d_R_mKFast;    ...
                d_L_hKFast;    d_R_hKFast;    ...
                d_L_syn_h;     d_R_syn_h;     ...
                d_L_syn_g;     d_R_syn_g];
    end



%% Initial Conditions 
rp = restPotential;
state0 = [ ...
    rp,                rp,                ...  % 1-2:   L/R vMembrane
    mNaFast_inf(rp),   mNaFast_inf(rp),   ...  % 3-4:   L/R mNaFast
    hNaFast_inf(rp),   hNaFast_inf(rp),   ...  % 5-6:   L/R hNaFast
    mNaPers_inf(rp),   mNaPers_inf(rp),   ...  % 7-8:   L/R mNaPers
    mCaFast_inf(rp),   mCaFast_inf(rp),   ...  % 9-10:  L/R mCaFast
    hCaFast_inf(rp),   hCaFast_inf(rp),   ...  % 11-12: L/R hCaFast
    mCaSlow_inf(rp),   mCaSlow_inf(rp),   ...  % 13-14: L/R mCaSlow
    hCaSlow_inf(rp),   hCaSlow_inf(rp),   ...  % 15-16: L/R hCaSlow
    mHypPol_inf(rp),   mHypPol_inf(rp),   ...  % 17-18: L/R mHypPol
    mKDelRect_inf(rp), mKDelRect_inf(rp), ...  % 19-20: L/R mKDelRect
    hKDelRect_inf(rp), hKDelRect_inf(rp), ...  % 21-22: L/R hKDelRect
    mKPers_inf(rp),    mKPers_inf(rp),    ...  % 23-24: L/R mKPers
    mKFast_inf(rp),    mKFast_inf(rp),    ...  % 25-26: L/R mKFast
    hKFast_inf(rp),    hKFast_inf(rp),    ...  % 27-28: L/R hKFast
    0, 0,                                 ...  % 29-30: L/R syn_h (zero at rest)
    0, 0                                  ...  % 31-32: L/R syn_g (zero at rest)
];

%% Solve
%tSpan = 0: timeStepSize: simulationTime;
options = odeset('RelTol',1e-6, 'AbsTol', 1e-8);

% Perform the numerical integration using ODE solver
%[timeStamps, solution] = ode15s(@odesys, tSpan, state0, options);
[timeStamps, solution] = ode15s(@odesys, [0, simulationTime], state0, options);
%[timeStamps, solution] = ode15s(@odesys, [0, simulationTime], state0);


%% Output
timeStamps  = timeStamps';
L_vMembrane = solution(:, 1)';
R_vMembrane = solution(:, 2)';

end