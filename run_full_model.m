%% Main file to run the full model
% (both within-host and human-vector systems)

tic
global P
set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 25 for 1X3, 20 for 1X2
%set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

RUN_constant = 1;
RUN_nonconstant = 0;

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
    % Generate the strategy via a combination of cubic splines
    % * Note that grid must match the splines grid *

    temp1 = importdata('basisMatrixNoKnots_1000_0.125.txt'); % choose from spline files
    CC1 = temp1.data(:,1);
    CC2 = temp1.data(:,2);
    CC3 = temp1.data(:,3);
    CC4 = temp1.data(:,4);

    % NB these weights calculated on [0,1000] with h = 0.5
    % beta = 12 optimal weights
    % w1 = -0.0947225580423;
    % w2 = 3.0996502357769;
    % w3 = -21.0071355096429;
    % w4 = 112.8119867596252;

    % beta = 14 optimal weights
    % w1 = 0.074769474519595;
    % w2 = 1.139206622742030;
    % w3 = -6.650902767084378;
    % w4 = 28.585817571349580;

    % % beta = 16 optimal weights, updated Nov 20
    w1 = 0.197629881402594;
    w2 = 0.168101173567905;
    w3 = -0.825428150237733;
    w4 = 2.193736391754480;

    % beta = 17 optimal weights
    % w1 = 0.274112296878810;
    % w2 = -0.037541099832533;
    % w3 = -0.017600836219813;
    % w4 = 0.060262196076839;

    % % beta = 16 with no immunity optimal weights
    % w1 = 0.199952113448875;
    % w2 = 0.196344830982411;
    % w3 = -0.818954166334971;
    % w4 = 1.970686551134877;

    CC = min(1,max(0,w1*CC1 + w2*CC2 + w3*CC3 + w4*CC4))';
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

%% Merozoite reproduction plotting
% P1 = h*trapz(I(1:length_infection_out,:).*repmat(gamma_fun(tau,h),1,length_infection_out)',2);
% P2 = h*trapz(I(1:length_infection_out,:),2)*P.mu;
% P3 = h*trapz(I(1:length_infection_out,:),2).*P.sigma.*(1-exp(-P.theta*A(1:length_infection_out)));
%
% R_M = (1-P.c).*P.beta*(P.p*B(1:length_infection_out,:)./(P.p*B(1:length_infection_out,:)+P.muM))...
%     .*(P1./(P1 + P2 + P3));

% figure(11);
% hold on;
% plot(x(1:length_infection_out)/24,R_M,'.-','LineWidth',3);
% xlabel('Time since infection (days)','Interpreter','latex');
% title('Effective Merozoite Number','Interpreter','latex');
% axis tight;

%%
toc