function b = withinhost_model_optimization(spline_weights)

%% numerical configuration
X_max = 280*24; % max time in days
tau_max = 20*24; % max 20 days
h = 0.25; % time/age step size in hours, same across all timescales

tau = (0:h:tau_max)';
ntau = length(tau);

% set model parameters via the baseline file (contains global variables)
baseline_parameter_set;

%% within-host model setup
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
%% Set the strategy based on the spline weights
temp1 = importdata('basisMatrixNoKnots_1000_0.5.txt'); % choose from spline files
CC1 = temp1.data(:,1);
CC2 = temp1.data(:,2);
CC3 = temp1.data(:,3);
CC4 = temp1.data(:,4);
%CC5 = temp1.data(:,5);
CC = min(1,max(0,spline_weights(1)*CC1 + spline_weights(2)*CC2 +...
    spline_weights(3)*CC3 + spline_weights(4)*CC4));
    %+ spline_weights(5)*CC5));

[~, ~, ~, ~, G, ~] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

b = -h*trapz(betaHV(G),1)/24; % return the cumulative infectiousness