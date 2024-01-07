%% Main file to run the full model
% (both within-host and human-vector systems)

tic
global P
set(0,'defaultaxesfontsize', 25);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none')
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')

%% numerical configuration
X_max = 300*24; % max time in days, max 200 days?
tau_max = 20*24; % max 20 days?
T_max = 200*24;
xV_max = 20*24;
h = 1; % time/age step size in hours, same across all timescales

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
% initially there are no merozoites or (developing/mature) gametocytes
B0 = P.Bstar; % scalar, nonzero
M0 = 0; % scalar, zero
I0 = ones(1,ntau); % I(0,tau), should be nonzero
I0(floor(48/h)+1:end) = 0; % I0 should be zero after 48 hours
initial_innoc = 0.06;
I0 = initial_innoc*I0/sum(I0);
% I0 uniform from zero to 48 hours approx.
IG0 = zeros(1,ntau); % IG(0,tau)
G0 = 0; % scalar, zero
A0 = 0; % scalar, zero

G_save = zeros(nx,131);
%% solve the within-host model for each value of P.c
for ii = 0:130
    P.c = 0.005*ii;
    [B, M, I, IG, G, A] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0);
    G_save(:,ii+1) = G;
end

%% Infectiousness plotting
figure;
imagesc(x/24,100*(0:0.01:0.65),betaHV(G_save)'); % Beta_HV(G(x)) heatmap
title('Infectiousness (\%)','Interpreter','latex');
xlabel('Time since infection (days)');
colormap jet;
colorbar;
clim([0 1]);
xlim([0 280]);
xticks([0 70 140 210 280]);
ytickformat('percentage');
set(gca,'YDir','normal');
%ylim([0 0.65]);
%%
toc