%% Main file to run the full model
% (both within-host and human-vector systems)

tic
global P
set(0,'defaultTextFontName', 'Arial');
set(0,'defaultaxesfontsize', 20); % 25 for 1X3, 20 for 1X2 figures
%set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% Choose constant or nonconstant investment strategy

RUN_constant = 0;
RUN_nonconstant = 1;
RUN_degree = 3; % degree of the polynomial spline (1,2,3,4 - 3 baseline)

%% numerical configuration
X_max = 1000*24; % max time in days
tau_max = 20*24; % max 20 days
T_max = 200*24;
xV_max = 20*24;
G_threshold = 1; % gametocyte threshold to end infection (doesn't impact dynamics)
h = 0.0625*2; % time/age step size in hours, same across all timescales

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

if RUN_constant
    % vector CC stores the strategy (proportion investmented in onward transmission)
    CC = P.c*ones(1,nx); % set the baseline constant investment strategy
else
    % Generate the strategy via a combination of the splines
    % * Note that grid must match the splines grid *
    if RUN_degree == 1
        temp1 = importdata('basisMatrixNoKnots_degree1_1000_0.125.txt'); % choose from spline files
        CC1 = temp1.data(:,1);
        CC2 = temp1.data(:,2);
        w1 = 0.094513570615327;
        w2 = -0.052554393112016;

        CC = min(1,max(0,w1*CC1 + w2*CC2))';
    elseif RUN_degree == 2
        temp1 = importdata('basisMatrixNoKnots_degree2_1000_0.125.txt'); % choose from spline files
        CC1 = temp1.data(:,1);
        CC2 = temp1.data(:,2);
        CC3 = temp1.data(:,3);
        w1 = 0.299711689528846;
        w2 = -0.374769344915881;
        w3 = 0.590730397030689;

        CC = min(1,max(0,w1*CC1 + w2*CC2 + w3*CC3))';
    elseif RUN_degree == 3
        %temp1 = importdata('basisMatrixKnot125_degree3_1000_0.125.txt'); %
        % knot options: 25, 75, 125
        temp1 = importdata('basisMatrixNoKnots_1000_0.125.txt');
        %standard
        CC1 = temp1.data(:,1);
        CC2 = temp1.data(:,2);
        CC3 = temp1.data(:,3);
        CC4 = temp1.data(:,4);

        % optimal weights with a knot at day 25, beta = 16
        % w1 = -23.842287422949948;
        % w2 = 0.484265317658271;
        % w3 = -0.455377758530129;
        % w4 = 0.400541070683154;

        % optimal weights with a knot at day 75, beta = 16
        % w1 = -0.222178862091422;
        % w2 = 0.512486123675708;
        % w3 = -0.513683492605781;
        % w4 = 0.478970226894363;

        % optimal weights with a knot at day 125, beta = 16
        % w1 = 0.002665299973821;
        % w2 = 0.390831987276833;
        % w3 = -0.376379753443364;
        % w4 = 0.371610178849970;

        % NB these weights are calculated on [0,1000] with h = 0.125
        % beta = 12 optimal weights
        % w1 = -0.0833658568039;
        % w2 = 3.7877052964783;
        % w3 = -34.5839214708051;
        % w4 = 249.3207851196597;

        % beta = 14 optimal weights
        % w1 = 0.072355007196402;
        % w2 = 1.212774197444973;
        % w3 = -7.336961107418777;
        % w4 = 32.746539009377429;

        % beta = 16 optimal weights -> STANDARD
        w1 = 0.197629881402594;
        w2 = 0.168101173567905;
        w3 = -0.825428150237733;
        w4 = 2.193736391754480;

        % beta = 17 optimal weights
        % w1 = 0.282994188967995;
        % w2 = -0.051947235680984;
        % w3 = -0.039442966387010;
        % w4 = 0.129188598123592;

        CC = min(1,max(0,w1*CC1 + w2*CC2 + w3*CC3 + w4*CC4))';
    elseif RUN_degree == 4
        %temp1 = importdata('basisMatrixKnot75_degree4_1000_0.125.txt');
        temp1 = importdata('basisMatrixNoKnots_degree4_1000_0.125.txt'); % choose from spline files
        CC1 = temp1.data(:,1);
        CC2 = temp1.data(:,2);
        CC3 = temp1.data(:,3);
        CC4 = temp1.data(:,4);
        CC5 = temp1.data(:,5);
        w1 = 0.220470836603252; % weights with no knots
        w2 = 0.296449499531031;
        w3 = -0.790096277376987;
        w4 = 1.050393368686625;
        w5 = 0.005878963403393;
        % w1 = -0.221087421864143; % weights with knot at day 75
        % w2 = 0.512198995086951;
        % w3 = -0.510947054140931;
        % w4 = 0.470198208196176;
        % w5 = 0.027510246532859;

        CC = min(1,max(0,w1*CC1 + w2*CC2 + w3*CC3 + w4*CC4 + w5*CC5))';
    else
        CC = P.c*ones(1,nx); % run constant strat if no valid degree
    end
end

% uncomment the following lines of code to perturb the given strategy
% sens_day = 580; % start time (in days) of the modification to strategy
% sens_length = 20; % number of days for which strategy is modified
% start_index = max(1,floor(sens_day*24/h));
% CC(1,start_index:(start_index+floor(sens_length*24/h)) )...
%     = 0.5*CC(1, start_index:(start_index+floor(sens_length*24/h)) );

[B, M, I, IG, G, A] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

%% solve the between-host/vector model
% HS0 = 100; % scalar
% HI0 = h*ones(1,nx); % HI(0,x), vector
% VS0 = 100; % scalar
% VI0 = zeros(1,nxV); % VI(t,xV) @ t = 0 (vector)
%
% [HS, HI, VS, VI] = human_vector_model(h, 0, X_max, T_max, xV_max, HS0, HI0, VS0, VI0, G);

%% Plotting and sim info

standard_plotting;
disp(['Cumulative infectiousness of the strategy: ',num2str(simps(0:h:X_max,betaHV(G))/24)]);

%%
toc