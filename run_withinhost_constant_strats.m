%% Compute the infectiousness profile for range of investment levels
tic

global P
set(0,'defaultaxesfontsize', 25);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none')
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesYGrid','on');

%% numerical configuration
X_max = 800*24; % max time in days (default around 700)
tau_max = 20*24; %  (default 20 days)
T_max = 200*24;
xV_max = 20*24;
h = 0.125; % time/age step size in hours, same across all timescales, 0.25 default
G_threshold = 1; % threshold for ending infection, 0.16341545 ~ 1% trans. prob.

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
initial_innoc = 0.06; % baseline: 0.06
I0 = initial_innoc*I0/(h*trapz(I0));
% I0 uniform from zero to 48 hours approx.
IG0 = zeros(1,ntau); % IG(0,tau)
G0 = 0; % scalar, zero
A0 = 0; % scalar, zero

invest_vec = 0:0.005:0.6; % vector of constant strategy percentages
%invest_vec = 0:0.0005:0.12; 
G_save = zeros(nx,length(invest_vec));
%% solve the within-host model for each value of P.c
for ii = 1:length(invest_vec)
    if mod(ii,50)==0
        disp([num2str(100*ii/length(invest_vec)),' % complete']);
    end
    P.c = invest_vec(ii);
    CC = P.c*ones(1,nx);
    [~, M, I, IG, G, A] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);
    G_save(:,ii) = G;
end

%% Infectiousness plotting
figure;
imagesc(x/24,100*invest_vec,betaHV(G_save)'); % Beta_HV(G(x)) heatmap
title('Infection Probability','Interpreter','latex');
xlabel('Time since infection (days)','Interpreter','latex');
colormap jet;
colorbar;
clim([0 1]);
set(gca,'ColorScale','linear')
xlim([0 650]);
%xticks([0 70 140 210 280 350]);
ytickformat('percentage');
set(gca,'YDir','normal');
ylabel('Transmission investment (\%)','Interpreter','latex');
ylim([0 60]);
grid off;
%% Gametocyte plotting heatmaps
% figure;
% imagesc(x/24,100*invest_vec,(G_save<G_threshold)');
% title('Presence of infection','Interpreter','latex');
% xlabel('Time since infection (days)','Interpreter','latex');
% colormap gray;
% xlim([0 650]);
% ytickformat('percentage');
% ylabel('Transmission investment (\%)','Interpreter','latex');
% set(gca,'YDir','normal');

% figure;
% imagesc(x/24,100*invest_vec,(G_save)');
% title('Gametocyte Count','Interpreter','latex');
% xlabel('Time since infection (days)','Interpreter','latex');
% colormap jetwhite;
% colorbar;
% xlim([0 500]);
% set(gca,'ColorScale','log')
% clim([1 10*10^4]);
% %xticks([0 70 140 210 280 350]);
% ytickformat('percentage');
% set(gca,'YDir','normal');
% ylabel('Transmission investment (\%)','Interpreter','latex');

%% Optimal strategy plotting
ac = floor(280*24/h)+1; % when we only want to integrate up to day ac
cum_inf1 = h*trapz(betaHV(G_save),1)/24;
%cum_inf1 = h*trapz(betaHV(G_save(1:ac,:)),1)/24;
figure(4);
hold on;
invest = 100*invest_vec;
plot(invest,cum_inf1,'-','Color',[0 0.4470 0.7410],'LineWidth',4);
psi = 1/105;
cum_inf2 = h*trapz(betaHV(G_save(1:ac,:)).*repmat(exp(-psi*x(1:ac)/24),1,length(invest_vec)),1)/24;
plot(invest,cum_inf2,':','Color',[0 0.4470 0.7410],'LineWidth',4);
% psi = 1/70;
% psi = 1/35;
xlim([0 max(invest)]);
ylim([0 400]);
xticks([0 10 20 30 40 50 60]);
xtickformat('percentage');
xtickangle(0);

% find the maxima for different values of psi (recovery rate)
[~, B] = max(cum_inf1);
scatter(invest(B),cum_inf1(B),200,'filled','k');
[~, B] = max(cum_inf2);
scatter(invest(B),cum_inf2(B),200,'filled','k');
% [~, B] = max(cum_inf3);
% scatter(invest(B),cum_inf3(B),200,'filled','k');
% [~, B] = max(cum_inf4);
% scatter(invest(B),cum_inf4(B),200,'filled','k');
ylabel('Cumulative Infectiousness ($f_1$)','Interpreter','latex');
xlabel('Transmission investment (\%)','Interpreter','latex');
%set(gca,'FontSize',35);
% legend('$\psi = 0$','$\psi = 1/105$','$\psi = 1/70$','$\psi = 1/35$',...
%     'Interpreter','latex','FontSize',35);
%% Optimal strategy versus time
% %ac = floor(350*24/h)+1;
% cum_inf1_time = h*cumtrapz(betaHV(G_save),1)/24;
% figure(6);
% hold on;
% invest = 100*invest_vec;
% %plot(x/24,cum_inf1_time(:,11),'LineWidth',3);
% plot(x/24,cum_inf1_time(:,B),'LineWidth',3);
% plot(x/24,cum_inf1_time(:,101),'LineWidth',3);
% plot(x/24,cum_inf1_time(:,181),'LineWidth',3);
% plot(x/24,cum_inf1_time(:,end),'LineWidth',3);
% 
% %colormap jet;
% %colorbar;
% xlim([0 X_max/24]);
% ylabel('Cumulative Infectiousness ($f_1$)','Interpreter','latex');
% xlabel('Time since infection (days)','Interpreter','latex');
% legend('$c = 4.4\%$','$c = 25\%$','$c = 45\%$','$c = 60\%$','Interpreter','latex','FontSize',25);

% % psi = 1/105;
% % int_range = (0:h:(ac-1)*h)/24;
% % cum_inf2 = h*sum(betaHV(G_save(1:ac,:)).*repmat(exp(-psi*int_range'),1,131),1)/24;
% % plot(invest,cum_inf2,'LineWidth',4);
% % psi = 1/70;
% % cum_inf3 = h*sum(betaHV(G_save(1:ac,:)).*repmat(exp(-psi*int_range'),1,131),1)/24;
% % plot(invest,cum_inf3,'LineWidth',4);
% % psi = 1/35;
% % cum_inf4 = h*sum(betaHV(G_save(1:ac,:)).*repmat(exp(-psi*int_range'),1,131),1)/24;
% % plot(invest,cum_inf4,'LineWidth',4);
% %xlim([0 max(invest)]);
% %xticks([0 10 20 30 40 50 60]);
% %xtickformat('percentage');
% %xtickangle(0);
%% Recovery time/Length of Infection plotting
% This plot is based on defining recovery as the last timet that there was
% > 1 gametocyte present in the host
temp_rec = zeros(1,length(invest));
for jj = 1:length(invest)
    rec_time = find(G_save(:,jj)>G_threshold,1,'last');
    if isempty(rec_time)
        temp_rec(jj) = X_max;
    else
        temp_rec(jj) = rec_time;
    end
end
figure(5);
hold on;
plot(invest,x(temp_rec)/24,'LineWidth',4);
%yline(280,'--','Color',[0 0.4470 0.7410],'LineWidth',4);
temp_duration = x(temp_rec)/24;
scatter(invest(B),temp_duration(B),200,'filled','k');
ylabel('Length of Infection (days)','Interpreter','latex');
xlabel('Transmission investment (\%)','Interpreter','latex');
xlim([0 max(invest)]);
ylim([0 650]);
%ylim([0 1.1*X_max/24]);
xticks([0 10 20 30 40 50 60]);
xtickformat('percentage');
xtickangle(0);
% legend('$\sigma = 0.75$','$\sigma = 0.7$','$\sigma = 0.65$','$\sigma = 0.6$',...
%     '$\sigma = 0.55$','$\sigma = 0.5$','Interpreter','latex','FontSize',25);
% legend('$\sigma \downarrow 20\%$','$\sigma (baseline)$','$\sigma \uparrow 10\%$','$\sigma \uparrow 20\%$',...
%     '$\sigma \uparrow 50\%$','$\sigma \uparrow 100\%$','Interpreter','latex','FontSize',25);
% legend('$\theta = 0.0000625$','$\theta = 0.000125$','$\theta = 0.0001875$','$\theta = 0.00025$',...
%     '$\theta = 0.0003125$','$\theta = 0.000375$','Interpreter','latex','FontSize',25);
%%
toc