%% Main file to run the full model 
% (both within-host and human-vector systems)

tic
global P
set(0,'defaultaxesfontsize', 25);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesYGrid','on');

%% numerical configuration
X_max = 800*24; % max time in days
tau_max = 20*24; % max 20 days?
T_max = 200*24;
xV_max = 20*24;
G_threshold = 1; % gametocyte threshold to end infection (doesn't impact dynamics)
h = 0.125; % time/age step size in hours, same across all timescales

x = (0:h:X_max)';
nx = length(x);
tau = (0:h:tau_max)';
ntau = length(tau);
t = (0:h:T_max)';
nt = length(t);
xV = (0:h:xV_max)';
nxV = length(xV);

% set model parameters via the baseline file (contains global variables)
baseline_parameter_set;

%% solve the within-host model
% initially there are no merozoites or (developing/mature) gametocytes
B0 = P.Bstar; % scalar, nonzero
M0 = 0; % scalar, zero
I0 = ones(1,ntau); % I(0,tau), should be nonzero
I0(floor(48/h)+1:end) = 0; % I0 should be zero after 48 hours
initial_innoc = 0.06; % baseline: 0.06
I0 = initial_innoc*I0/(h*trapz(I0));
% I0 uniform from zero to 48 hours approx.
IG0 = zeros(1,ntau); % IG(0,tau)
G0 = 0; % scalar, zero
A0 = 0; % scalar, zero

% NB: ordering of independent variables is I(x,tau), IG(x,tau)


%% Set the parasite investment strategy and run the within-host model

% vector CC stores the strategy (proportion investmented in onward transmission)
CC = P.c*ones(1,nx); % set the baseline constant investment strategy

% Generate the strategy via a combination of cubic splines
% temp1 = importdata('basisMatrixNoKnots.txt'); % choose from spline files
% CC1 = temp1.data(:,1);
% CC2 = temp1.data(:,2);
% CC3 = temp1.data(:,3);
% CC4 = temp1.data(:,4);
% w1 = 0.199944393818933; 
% w2 = 0.197075562766502; 
% w3 = -0.453310918950100; 
% w4 = 0.714969288272976; % weights for basis splines
% CC = max(0,w1*CC1 + w2*CC2 + w3*CC3 + w4*CC4); 


% uncomment the following lines of code to perturb the constant strategy
% sens_day = 3; % start time (in days) of the modification to strategy
% sens_length = 2; % number of days for which strategy is modified
% start_index = max(1,floor(sens_day*24/h));
% CC(1,start_index:(start_index+floor(sens_length*24/h)) )...
%     = 2*CC(1, start_index:(start_index+floor(sens_length*24/h)) );

[B, M, I, IG, G, A] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

%% solve the between-host/vector model
% HS0 = 100; % scalar
% HI0 = h*ones(1,nx); % HI(0,x), vector
% VS0 = 100; % scalar
% VI0 = zeros(1,nxV); % VI(t,xV) @ t = 0 (vector)
% 
% [HS, HI, VS, VI] = human_vector_model(h, 0, X_max, T_max, xV_max, HS0, HI0, VS0, VI0, G);

%%
standard_plotting;

%%
figure(7);
hold on;
title('$I(x,0)$','Interpreter','latex');
plot(x/24,I(:,1),lineStyle,'LineWidth',3); 
xlim([0 10]);

%%
figure(8);
hold on;
title('$IG(x,0)$','Interpreter','latex');
plot(x/24,IG(:,1),lineStyle,'LineWidth',3); 
xlim([0 10]);

%%
toc