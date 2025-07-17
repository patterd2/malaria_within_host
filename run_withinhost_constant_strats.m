%% Compute the infectiousness profile for range of investment levels
tic

global P
set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% Plotting toggles
HEATMAP = 0; % Toggle for heatmap plotting, slow at high resolutions

%% numerical configuration
X_max = 750*24; % max time in days
tau_max = 20*24; %  (default 20 days)
T_max = 200*24;
xV_max = 20*24;
h = 0.0625*2; % time/age step size in hours, same across all timescales, 0.25 default
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

invest_vec = 0.01:0.01:0.6; %  vector of constant strategy percentages
%invest_vec = 0.005:0.005:0.02; 0.0205:0.005:0.1;
G_save = zeros(nx,length(invest_vec));
%% solve the within-host model for each value of P.c
for ii = 1:length(invest_vec)
    if mod(ii,20)==0
        disp([num2str(100*ii/length(invest_vec)),' % complete']);
    end
    P.c = invest_vec(ii);
    CC = P.c*ones(1,nx);
    [~, M, I, IG, G, A] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);
    G_save(:,ii) = G;
end

%% Host infectiousness heatmap plotting
if HEATMAP
    figure(3);
    imagesc(x/24,100*invest_vec,betaHV(G_save)'); % Beta_HV(G(x)) heatmap
    LimitsX = xlim; LimitsY = ylim;
    title('A. host infectiousness','FontWeight','normal',...
        'HorizontalAlignment','left','position', [LimitsX(1), LimitsY(2)]);
    xlabel('infection age x (days)');
    colormap jet;
    colorbar;
    clim([0 1]);
    set(gca,'ColorScale','linear')
    xlim([0 650]);
    ytickformat('percentage');
    set(gca,'YDir','normal');
    ylabel('transmission investment (c)');
    ylim([0.5 60]);
    grid off;
    set(gca,'TickDir','out');
    yline(4.4,'LineWidth',4,'Color',[0.5 0.5 0.5],'Alpha',1);
    legend('\color{gray} optimal strategy','Interpreter','tex');
    legend('boxoff')
    box off;
end
%% Recovery time/Length of Infection plotting
% calculate the last time the gametocyte population is above the threshold
invest = 100*invest_vec;
temp_rec = zeros(1,length(invest));
G_tmp = zeros(size(G_save));
for jj = 1:length(invest)
    rec_time = find(G_save(:,jj) > G_threshold,1,'last');
    if isempty(rec_time)
        if G_save(:,jj) > G_threshold
            temp_rec(jj) = X_max;
        else
            temp_rec(jj) = 1;
        end
    else
        % record the gametocytes up to the threshold for fitness calc
        G_tmp(1:rec_time,jj) = G_save(1:rec_time,jj);
        temp_rec(jj) = rec_time;
    end
end

%% Optimal strategy plotting
ac = floor(280*24/h)+1; % when we only want to integrate up to day ac, i.e. no immunity case
if P.sigma == 0
    cum_inf1 = trapz(x(1:ac),betaHV(G_save(1:ac,:)))/24;
else
    cum_inf1 = h*trapz(betaHV(G_tmp))/24; % fitness calculation
end

%% Plot fitness versus transmission investment parameter c
figure(4);
hold on;
plot(invest,cum_inf1,'-','Color',[0 0.4470 0.7410],'LineWidth',4);
psi = 1/105;
cum_inf2 = simps(x(1:ac),betaHV(G_save(1:ac,:)).*repmat(exp(-psi*x(1:ac)/24),1,length(invest_vec)),1)/24;
plot(invest,cum_inf2,':','Color',[0 0.4470 0.7410],'LineWidth',4);
% psi = 1/70; 
% psi = 1/35;
xlim([0.0 max(invest)]);
ylim([0 300]);
xticks([0 20 40 60]);
xlim([0 60]);
%yticks([0 100 200 300]);
xtickformat('percentage');
xtickangle(0);

% find the maxima for different values of psi (recovery rate)
[~, B1] = max(cum_inf1);
scatter(invest(B1),cum_inf1(B1),200,'filled','k');
[~, B2] = max(cum_inf2);
scatter(invest(B2),cum_inf2(B2),200,'filled','k');
%[~, B] = max(cum_inf3);
%scatter(invest(B),cum_inf3(B),200,'filled','k');
%[~, B] = max(cum_inf4);
%scatter(invest(B),cum_inf4(B),200,'filled','k');
ylabel('cumulative infectiousness');
xlabel('costant transmission investment (c)');
set(gca,'TickDir','out');
box off;
LimitsX = xlim; LimitsY = ylim;
title('B. fitness','FontWeight','normal',...
    'HorizontalAlignment','left','position', [LimitsX(1), LimitsY(2)]);

%% Length of infection as a function of transmission investment plot
figure(5);
hold on;
plot(invest,x(temp_rec)/24,'--','Color',[0.47,0.67,0.19],'LineWidth',4);
%yline(280,'--','Color',[0 0.4470 0.7410],'LineWidth',4,'Alpha',1);
%yline(1/psi,':','Color',[0 0.4470 0.7410],'LineWidth',4,'Alpha',1);
temp_duration = x(temp_rec)/24;
scatter(invest(B1),temp_duration(B1),200,'filled','k');
ylabel('infection duration (days)');
xlabel('constant transmission investment (c)');
xlim([0.0 max(invest)]);
ylim([0 750]);
%ylim([0 1.1*X_max/24]);
xticks([0 20 40 60]);
xlim([0 60]);
xtickformat('percentage');
xtickangle(0);
box off;
set(gca,'TickDir','out');
%legend('immune feedback','no immunity','exp.recovery rate');
%legend('boxoff')
LimitsX = xlim; LimitsY = ylim;
title('C. duration of infection','FontWeight','normal',...
    'HorizontalAlignment','left','position', [LimitsX(1), LimitsY(2)]);
%%
toc