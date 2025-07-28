function b = withinhost_model_optimization(spline_weights)


%% numerical configuration
X_max = 1000*24; % max time in days
tau_max = 20*24; % max 20 days
h = 0.0625*2; % time/age step size in hours, same across all timescales

tau = (0:h:tau_max)';
ntau = length(tau);
x = (0:h:X_max)';
nx = length(x);
ac = floor(280*24/h)+1; % for no immune fitness calculations
psi = 1/105; % constant recovery rate in no immunity case

G_threshold = 1; % gametocyte threshold to end infection for fitness calc

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
if isscalar(spline_weights)
    CC = spline_weights*ones(1,nx);
elseif length(spline_weights) == 2
    temp1 = importdata('basisMatrixNoKnots_degree1_1000_0.125.txt'); % choose linear splines
    CC1 = temp1.data(:,1);
    CC2 = temp1.data(:,2);
    CC = min(1,max(0,spline_weights(1)*CC1 + spline_weights(2)*CC2));
elseif length(spline_weights) == 3
    temp1 = importdata('basisMatrixNoKnots_degree2_1000_0.125.txt'); % choose linear splines
    CC1 = temp1.data(:,1);
    CC2 = temp1.data(:,2);
    CC3 = temp1.data(:,3);
    CC = min(1,max(0,spline_weights(1)*CC1 + spline_weights(2)*CC2 + spline_weights(3)*CC3));
elseif length(spline_weights) == 4
    %temp1 = importdata('basisMatrixKnot125_degree3_1000_0.125.txt');
    %temp1 = importdata('basisMatrixKnot25_degree3_1000_0.125.txt');
    temp1 = importdata('basisMatrixKnot75_degree3_1000_0.125.txt');
    %temp1 = importdata('basisMatrixNoKnots_1000_0.125.txt'); % no knot
    CC1 = temp1.data(:,1);
    CC2 = temp1.data(:,2);
    CC3 = temp1.data(:,3);
    CC4 = temp1.data(:,4);
    CC = min(1,max(0,spline_weights(1)*CC1 + spline_weights(2)*CC2 +...
        spline_weights(3)*CC3 + spline_weights(4)*CC4));
elseif length(spline_weights) == 5
    temp1 = importdata('basisMatrixKnot75_degree4_1000_0.125.txt'); % choose from spline files
    %temp1 = importdata('basisMatrixNoKnots_degree4_1000_0.125.txt'); % choose linear splines
    CC1 = temp1.data(:,1);
    CC2 = temp1.data(:,2);
    CC3 = temp1.data(:,3);
    CC4 = temp1.data(:,4);
    CC5 = temp1.data(:,5);
    CC = min(1,max(0,spline_weights(1)*CC1 + spline_weights(2)*CC2 +...
        spline_weights(3)*CC3 + spline_weights(4)*CC4 + spline_weights(5)*CC5));
end

[~, ~, ~, ~, G, ~] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

length_infection = find(G>G_threshold,1,'last');
if isempty(length_infection)
    length_infection = length(x);
end
if P.sigma == 0
    b = -simps(x(1:ac),betaHV(G(1:ac)).*exp(-psi*x(1:ac)/24))/24;
    % return the cumulative infectiousness in no immunity case, f_prev
else
    b = -simps(x(1:length_infection),betaHV(G(1:length_infection)))/24; 
    % return the cumulative infectiousness/fitness, f
end