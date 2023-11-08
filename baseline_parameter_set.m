global P
% values from Greischar et al., need to check units
%% within-host model parameters
P.c = 0.2; % parasite investment fraction
P.lambda = (2*10^(5))/24;
P.Bstar = 5*10^6;
P.K = P.lambda*(P.Bstar)/(P.lambda - P.mu*P.Bstar); % ?
P.p = (8.35*10^(-6))/24;
P.beta = 16;
P.mu = (1/120)/24;
P.muM = 200/24; 
P.muG = 0.5/24;

% immune activation function parameters
P.sigma = 0; % sigma = 0 turns off immune system
P.IT = 0.25; % immune activation sigmoid threshold (phi)
P.s = 1; % immune activation sigmoid slope (phi)
P.theta = 1; 
% need to tune these parameters...

%% Human and vector parameters
P.b = 0.35;
P.N = 2000;
P.betaVH = 0.05;
P.deltaA = 0.093;