global P
% values from Greischar et al., need to check units
%% within-host model parameters
P.c = 0.044; % parasite investment fraction
P.mu = (1/120)/24;
P.lambda = (2*10^(5))/24;
P.Bstar = 5*10^6;
P.K = P.lambda*(P.Bstar)/(P.lambda - P.mu*P.Bstar); % ?
P.p = (8.35*10^(-6))/24; % (8.35*10^(-6))/24 baseline
P.beta = 16; % 16
P.muM = 200/24;
P.muG = 0.5/24;
P.Bbar = P.lambda/(P.mu + P.lambda/P.K);
P.R0 = ((P.beta*(1-P.c))*P.p/(P.muM + P.p*P.Bbar));

% immune activation function parameters
P.sigma = 0.55/24; % sigma = 0 turns off immune system, 0.55/24 baseline
P.IT = 2; % immune activation sigmoid threshold (phi)
P.s = 1; % immune activation sigmoid slope (phi), Heaviside at zero
P.theta = 0.00025; % set very high to have constant immune removal
P.muA = 0.0/24; % set to zero to turn off immune relaxation, need to "tune"

%% Human and vector parameters
P.b = 0.35/24;
P.N = 2000;
P.betaVH = 0.05;
P.deltaA = 0.093/24;