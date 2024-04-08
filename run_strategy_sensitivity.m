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

%% numerical configuration
X_max = 700*24; % max time in days, max 300 days?
tau_max = 20*24; % max 20 days?
T_max = 200*24;
xV_max = 20*24;
h = 0.5; % time/age step size in hours, same across all timescales

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
CC = P.c*ones(1,nx); % set the investment strategy, 10% is optimal
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
    CC = P.c*ones(1,nx); % set the investment strategy
    CC(1,1+(ii-1)*2*24/h:ii*2*24/h) = 0.5*CC(1,1+(ii-1)*2*24/h:ii*2*24/h);
    [~, ~, ~, ~, G, ~] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);
    cum_inf1_dec(ii) = h*sum(betaHV(G),1)/24;
    length_infection_dec(ii) = find(G>1,1,'last');
    CC = P.c*ones(1,nx); % set the investment strategy
    CC(1,1+(ii-1)*2*24/h:ii*2*24/h) = 1.5*CC(1,1+(ii-1)*2*24/h:ii*2*24/h);
    [~, ~, ~, ~, G, ~] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);
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