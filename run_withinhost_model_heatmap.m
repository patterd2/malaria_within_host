%% Compute the infectiousness profile for range of investment levels
tic

global P
set(0,'defaultaxesfontsize', 25);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none')
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')

%% numerical configuration
X_max = 700*24; % max time in days, max 200 days?
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
I0 = initial_innoc*I0/(h*trapz(I0));
% I0 uniform from zero to 48 hours approx.
IG0 = zeros(1,ntau); % IG(0,tau)
G0 = 0; % scalar, zero
A0 = 0; % scalar, zero

G_save = zeros(nx,131);
%% solve the within-host model for each value of P.c
for ii = 0:130
    P.c = 0.005*ii;
    [~, M, I, IG, G, A] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0);
    G_save(:,ii+1) = G;
end

%% Infectiousness plotting
figure;
imagesc(x/24,100*(0:0.005:0.65),betaHV(G_save)'); % Beta_HV(G(x)) heatmap
title('Infectiousness (\%)','Interpreter','latex');
xlabel('Time since infection (days)');
colormap jet;
colorbar;
clim([0 1]);
xlim([0 350]);
xticks([0 70 140 210 280 350]);
ytickformat('percentage');
set(gca,'YDir','normal');
%ylim([0 0.65]);
%% Gametocyte plotting heatmaps
figure;
imagesc(x/24,100*(0:0.005:0.65),(G_save<1)');
title('Gametocytes above threshold?','Interpreter','latex');
xlabel('Time since infection (days)');
colormap gray;
%colorbar;
xlim([0 770]);
%xticks([0 70 140 210 280 350]);
ytickformat('percentage');
ylabel('Transmission investment (\%)','Interpreter','latex');
set(gca,'YDir','normal');

figure;
imagesc(x/24,100*(0:0.005:0.65),(G_save)');
title('Gametocyte Count','Interpreter','latex');
xlabel('Time since infection (days)','Interpreter','latex');
colormap jet;
colorbar;
xlim([0 210]);
%xticks([0 70 140 210 280 350]);
ytickformat('percentage');
set(gca,'YDir','normal');
ylabel('Transmission investment (\%)','Interpreter','latex');

%% Optimal strategy plotting
ac = floor(350*24/h)+1;
cum_inf1 = h*sum(betaHV(G_save(1:ac,:)),1)/24;
figure;
invest = 100*(0:0.005:0.65);
plot(invest,cum_inf1,'LineWidth',4);
hold on;
% psi = 1/105;
% int_range = (0:h:(ac-1)*h)/24;
% cum_inf2 = h*sum(betaHV(G_save(1:ac,:)).*repmat(exp(-psi*int_range'),1,131),1)/24;
% plot(invest,cum_inf2,'LineWidth',4);
% psi = 1/70;
% cum_inf3 = h*sum(betaHV(G_save(1:ac,:)).*repmat(exp(-psi*int_range'),1,131),1)/24;
% plot(invest,cum_inf3,'LineWidth',4);
% psi = 1/35;
% cum_inf4 = h*sum(betaHV(G_save(1:ac,:)).*repmat(exp(-psi*int_range'),1,131),1)/24;
% plot(invest,cum_inf4,'LineWidth',4);
xlim([0 65]);
xticks([0 10 20 30 40 50 60]);
xtickformat('percentage');
% find the maxima for different values of psi (recovery rate)
[~, B] = max(cum_inf1);
scatter(invest(B),cum_inf1(B),200,'filled','k');
% [~, B] = max(cum_inf2);
% scatter(invest(B),cum_inf2(B),200,'filled','k');
% [~, B] = max(cum_inf3);
% scatter(invest(B),cum_inf3(B),200,'filled','k');
% [~, B] = max(cum_inf4);
% scatter(invest(B),cum_inf4(B),200,'filled','k');
ylabel('cumulative infectiousness ($f_1$)','Interpreter','latex');
xlabel('Transmission investment (\%)','Interpreter','latex');
set(gca,'FontSize',35);

% legend('$\psi = 0$','$\psi = 1/105$','$\psi = 1/70$','$\psi = 1/35$',...
%     'Interpreter','latex','FontSize',35);
%% Recovery time/Length of Infection plotting
% This plot is based on defining recovery as the last timet that there was
% > 1 gametocyte present in the host
temp_rec = zeros(1,length(invest));
for jj = 1:length(invest)
    rec_time = find(G_save(:,jj)>1,1,'last');
    if isempty(rec_time)
        temp_rec(jj) = X_max;
    else
        temp_rec(jj) = rec_time;
    end
end
figure;
plot(invest,temp_rec/24,'LineWidth',3);
ylabel('Length of Infection (days)','Interpreter','latex');
xlabel('Transmission investment (\%)','Interpreter','latex');
xlim([0 max(invest)]);
ylim([0 X_max*1.1/24]);
%%
toc