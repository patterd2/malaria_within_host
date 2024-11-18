%% Bang-bang-bang sensitivity analysis of strategies
% (both within-host and human-vector systems
tic
global P
set(0,'defaultaxesfontsize', 25);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesYGrid','on');

RUN_constant = 1;

%% numerical configuration
X_max = 700*24; % max time in days, max 300 days?
tau_max = 20*24; % max 20 days?
T_max = 200*24;
xV_max = 20*24;
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

%% Calculate cumulative infectiousness baseline
if RUN_constant
    % vector CC stores the strategy (proportion investmented in onward transmission)
    CC = P.c*ones(1,nx); % set the baseline constant investment strategy
else
    % Generate the strategy via a combination of cubic splines
    % * Note that grid must match the splines grid *
    
    temp1 = importdata('basisMatrixNoKnots_1000_0.5.txt'); % choose from spline files
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

    % % beta = 16 optimal weights
    w1 = 0.199952113448875;
    w2 = 0.196344830982411;
    w3 = -0.818954166334971;
    w4 = 1.970686551134877;

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
[~, ~, ~, ~, G, ~] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);
cum_inf1_baseline = h*trapz(betaHV(G),1)/24;
length_infection_baseline = find(G>1,1,'last');
length_infection_baseline = x(length_infection_baseline)/24; % infection length in days
%%
sim_days = (X_max/48);
cum_inf1_inc = zeros(1,sim_days);
cum_inf1_dec = zeros(1,sim_days);
length_infection_inc = zeros(1,sim_days);
length_infection_dec = zeros(1,sim_days);
for ii = 1:sim_days % changing strategy every two days
    CC_dec = CC;
    CC_dec(1,1+(ii-1)*2*24/h:ii*2*24/h) = 0.5*CC(1,1+(ii-1)*2*24/h:ii*2*24/h);
    [~, ~, ~, ~, G, ~] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC_dec);
    cum_inf1_dec(ii) = h*sum(betaHV(G),1)/24;
    length_infection_dec(ii) = find(G>1,1,'last');
    CC_inc = CC; 
    CC_inc(1,1+(ii-1)*2*24/h:ii*2*24/h) = 2.0*CC(1,1+(ii-1)*2*24/h:ii*2*24/h);
    [~, ~, ~, ~, G, ~] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC_inc);
    cum_inf1_inc(ii) = h*sum(betaHV(G),1)/24;
    length_infection_inc(ii) = find(G>1,1,'last');
end
%%
figure;
plot(2*(1:sim_days),cum_inf1_dec-cum_inf1_baseline,'LineWidth',3);
hold on;
plot(2*(1:sim_days),cum_inf1_inc-cum_inf1_baseline,'LineWidth',3);
%plot(2*(1:sim_days),cum_inf1_baseline*ones(1,length(1:sim_days)),':','LineWidth',3);
xlim([0 700]);
hold on;
xlabel('Day','Interpreter','latex');
ylabel('Cumulative Infectiousness ($f_1$)','Interpreter','latex');
legend('50\% decrease','50\% increase');
%%
figure;
plot(2*(1:sim_days),x(length_infection_dec)/24 - length_infection_baseline,'LineWidth',3);
hold on;
plot(2*(1:sim_days),x(length_infection_inc)/24 - length_infection_baseline,'LineWidth',3);
%plot(2*(1:sim_days),length_infection_baseline*ones(1,length(1:sim_days)),'--','LineWidth',3);
xlim([0 700]);
xlabel('Day','Interpreter','latex');
ylabel('Length of infection (days)','Interpreter','latex');
legend('50\% decrease','50\% increase','Location','northwest');
%%
% figure;
% plot(h*(1:sim_days),length_infection_inc,'LineWidth',3);
% title('Sensitivity to 50% increases');
% xlabel('Day','Interpreter','latex');
% ylabel('Length of infection',Interpreter','latex');
%%
toc