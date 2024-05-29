%% numerical configuration
X_max = 600*24; % max time in days, max 300 days?
tau_max = 20*24; % max 20 days?
T_max = 200*24;
xV_max = 20*24;
h = 0.03125; % time/age step size in hours, same across all timescales

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
CC = P.c*ones(1,nx); % set the investment strategy

%% Calculate cumulative infectiousness baseline
tic
[~, ~, ~, ~, G, ~] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);
toc