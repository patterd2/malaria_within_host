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
X_max = 700*24; % max time in days, max 300 days?
tau_max = 20*24; % max 20 days?
T_max = 200*24;
xV_max = 20*24;
G_threshold = 1; % gametocyte threshold to end infection
h = 0.25; % time/age step size in hours, same across all timescales

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
initial_innoc = 0.06;
I0 = initial_innoc*I0/(h*trapz(I0));
% I0 uniform from zero to 48 hours approx.
IG0 = zeros(1,ntau); % IG(0,tau)
G0 = 0; % scalar, zero
A0 = 0; % scalar, zero

% NB: ordering of independent variables is I(x,tau), IG(x,tau)

CC = P.c*ones(1,nx); % set the baseline constant investment strategy

sens_day = 2; % start day of the modification to strategy
sens_length = 3; % number of days for which strategy is modified
CC(1,1+(sens_day-1)*sens_length*24/h:sens_day*sens_length*24/h)...
    = 2*CC(1, 1+(sens_day-1)*sens_length*24/h:sens_day*sens_length*24/h );
% sensitivity checks: increase constant strat by 50% in chosen 48 hours

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