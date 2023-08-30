global P
% values from Greischar et al., need to check units
%% within-host model parameters
P.c = 0.5; % parasite investment fraction
P.lambda = 1;%(2*10^(5));
P.K = 100; % ?
P.p = 1;%(8.35*10^(-6));
P.beta = 16;
P.mu = 1/120;
P.muM = 200; 
P.muG = 0.5;
% immune activation function parameters
P.sigma = 0; % sigma = 0 turns off immune system
P.IT = 1;
P.s = 1;
P.theta = 1; % need to tune these parameters...

%% Human and vector parameters
P.b = 0.35;
P.N = 2000;
P.betaVH = 0.05;
P.deltaA = 0.093;