global P
% values from Greischar et al., need to check units
%% within-host model parameters
P.c = 0.044; % parasite investment fraction (baseline parameters 0.044)
P.mu = (1/120)/24; % 1/120 baseline (daily rate)
P.lambda = (2*10^(5))/24; % 2*10^(5) baseline (daily rate)
P.Bstar = 5*10^6; % 5*10^6 baseline
P.K = P.lambda*(P.Bstar)/(P.lambda - P.mu*P.Bstar);
P.p = (8.35*10^(-6))/24; % (8.35*10^(-6))/24 baseline
P.beta = 16; % 16 baseline
P.muM = 200/24; % 200 baseline
P.muG = 0.5/24; % 0.5 baseline
P.Bbar = P.lambda/(P.mu + P.lambda/P.K);

% immune activation function parameters
P.sigma = 0.55/24; % sigma = 0 turns off immune system, 0.55/24 baseline
P.IT = 2; % immune activation sigmoid threshold (phi), 2 baseline
P.s = 1; % immune activation sigmoid slope (phi), Heaviside at zero, baseline 1
P.theta = 0.00025; % 0.00025 baseline
P.muA = 0.0/24; % set to zero to turn off immune relaxation, 0 baseline

%% Human and vector parameters
P.b = 0.35/24;
P.N = 2000;
P.betaVH = 0.05;
P.deltaA = 0.093/24;