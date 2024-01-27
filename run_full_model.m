%% Main file to run the full model 
% (both within-host and human-vector systems)

tic
global P
set(0,'defaultaxesfontsize', 25);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesYGrid','on');

%% numerical configuration
X_max = 300*24; % max time in days, max 200 days?
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

%% solve the within-host model
% initially there are no merozoites or (developing/mature) gametocytes
B0 = P.Bstar; % scalar, nonzero
M0 = 0; % scalar, zero
I0 = ones(1,ntau); % I(0,tau), should be nonzero
I0(floor(48/h)+1:end) = 0; % I0 should be zero after 48 hours
initial_innoc = 0.06;
I0 = initial_innoc*I0/sum(I0);
% I0 uniform from zero to 48 hours approx.
IG0 = zeros(1,ntau); % IG(0,tau)
G0 = 0; % scalar, zero
A0 = 0; % scalar, zero

% NB: ordering of independent variables is I(x,tau), IG(x,tau)

[B, M, I, IG, G, A] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0);
%% solve the between-host/vector model
HS0 = 100; % scalar
HI0 = h*ones(1,nx); % HI(0,x), vector
VS0 = 100; % scalar
VI0 = zeros(1,nxV); % VI(t,xV) @ t = 0 (vector)

[HS, HI, VS, VI] = human_vector_model(h, 0, X_max, T_max, xV_max, HS0, HI0, VS0, VI0, G);
%% plot within-host dynamics
figure(1);
subplot(1,3,1), plot(x/24,B,'--','LineWidth',3); hold on;
title('$B(x)$','Interpreter','latex');
hold on;
subplot(1,3,2), plot(x/24,M,'--','LineWidth',3); hold on;
title('$M(x)$','Interpreter','latex');
xlabel('Age of infection (x) [days]');
subplot(1,3,3), plot(x/24,G,'--','LineWidth',3); hold on;
title('$G(x)$','Interpreter','latex');
%legend('$B(x)$ uninfected red blood cells','$M(x)$ merozoites','$G(x)$ gametocytes');

%%
% figure;
% plot(x/24,I(:,1)); % I(x,tau)
% hold on;
% plot(x/24,I(:,floor(length(x)/6)));
% plot(x/24,I(:,floor(length(x)/3)));
% plot(x/24,I(:,floor(2*length(x)/3)));
% plot(x/24,I(:,floor(length(x))));
% title('Infection dynamics (asexual stage): $I(x,\tau)$','Interpreter','latex');
% xlabel('Time since infection (days)');
% legend('$I(x,0)$','$I(x,0.5)$','$I(x,1)$','$I(x,2)$','$I(x,3)$');
%axis tight;

%% 
% figure(2);
% plot(x/24,h*sum(I,2),'LineWidth',3); % I(x,tau)
% hold on;
% title('Infection dynamics (asexual stage): $\int I(x,\tau) \, d\tau$','Interpreter','latex');
% xlabel('Time since infection (days)');
% axis tight;
% legend('$h = 2$','$h = 1$','$h = 0.5$','$h = 0.25$','FontSize',35);
% legend('$c = 0.05$','$c = 0.4$','FontSize',35);

%% Infectiousness plotting
% figure(3);
% plot(x/24,betaHV(G),'LineWidth',3); % Beta_HV(G(x))
% hold on;
% title('Infectiousness (\%)','Interpreter','latex');
% xlabel('Time since infection (days)');
% ylim([0 1]);

%% Combined plotting
figure(4);
subplot(2,2,1), plot(x/24,h*sum(I,2),'--','LineWidth',3); hold on;
axis tight;
title('$\int I(x,\tau) \, d\tau$ (asexual stage)','Interpreter','latex');
hold on;
subplot(2,2,2), plot(x/24,G,'--','LineWidth',3); hold on;
title('$G(x)$','Interpreter','latex');
subplot(2,2,3), plot(x/24,A,'--','LineWidth',3); hold on;
xlabel('Age of infection (x) [days]');
title('$A(x)$','Interpreter','latex');
subplot(2,2,4), plot(x/24,betaHV(G),'--','LineWidth',3); hold on;
title('Infectiousness (\%)','Interpreter','latex');
xlabel('Age of infection (x) [days]');
% legend('$c = 0.05$','$c = 0.4$','FontSize',35);
%legend('$\mu_A = 0$','$\mu_A = 0.01/24$','$\mu_A = 0.1/24$','FontSize',35);
%% Immune dynamics plotting
%figure(5);
% subplot(2,1,1), plot(x/24,phi(h*sum(I,2),P.IT,P.s),'LineWidth',3); hold on;
% axis tight;
% xlabel('Age of infection (x) [days]');
% title('$\phi(\int I(x,\tau) \, d\tau)$','Interpreter','latex');
% hold on;
% subplot(2,1,2), 
%plot(x/24,A,'LineWidth',3); hold on;
%title('$A(x)$','Interpreter','latex');
%xlabel('Age of infection (x) [days]');
%legend('$\sigma = 0$','$\sigma = 0.1$',...
%    '$\sigma = 0.25$','$\sigma = 0.5$','FontSize',35);
%% Host infection plotting
% figure;
% plot(t/24,sum(HI,2),'LineWidth',3);
% title('$\int H_I(t,x) \, dx$','Interpreter','latex');
% xlabel('Time (days)');
%%

% figure;
% plot(x/24,IG(:,1)); % I(x,tau)
% hold on;
% plot(x/24,IG(:,floor(length(x)/6)));
% plot(x/24,IG(:,floor(length(x)/3)));
% plot(x/24,IG(:,floor(2*length(x)/3)));
% plot(x/24,IG(:,floor(length(x))));
% title('Infection dynamics (sexual stage): $I_G(x,\tau)$');
% xlabel('Time since infection (days)');
% legend('$I_G(x,0)$','$I_G(x,0.5)$','$I_G(x,1)$','$I_G(x,2)$','$I_G(x,3)$');
% axis tight;
% grid on;

%%
% figure;
% plot(x/24,A,'LineWidth',3);
% hold on
% plot(x/24,P.sigma*(1-exp(-P.theta*A)),'LineWidth',3);
% hold on;
% xlabel('Time since infection (days)');
% legend('$A(x)$ immune activation level','$\sigma(1 - \exp(\theta A(x)))$ infection removal rate');
% grid on;

%%
% figure;
% imagesc(x/24,x/24,I);
% title('$I(x,\tau)$');
% set(gca,'YDir','normal');
% colorbar;
% clim([min(min(I)) max(max(I))+1]);
% xlabel('$x$');
% ylabel('$\tau$');
% grid on;

% figure;
% imagesc(x/24,x/24,IG);
% title('$I_G(x,\tau_G)$');
% set(gca,'YDir','normal');
% colorbar;
% clim([min(min(IG)) max(max(IG))+1]);
% xlabel('$x$');
% ylabel('$\tau_G$');
% grid on;

%% plot human-vector dynamics
% figure_setups;
% plot(x/24,HS,'LineWidth',2);
% hold on;
% plot(x/24,VS,'LineWidth',2);
% xlabel('Time (days)');
% legend('$H_S(t)$ susceptible humans','$V_S(t)$ susceptible vectors');
% grid on;
% 
% figure;
% imagesc(x/24,x/24,HI);
% title('$H_I(t,x)$');
% set(gca,'YDir','normal');
% colorbar;
% clim([min(min(HI)) max(max(HI))+1]);
% xlabel('$t$');
% ylabel('$x$');
% grid on;
% 
% figure;
% imagesc(x/24,x/24,VI);
% title('$V_I(t,\tau_V)$');
% set(gca,'YDir','normal');
% colorbar;
% clim([min(min(VI)) max(max(VI))+1]);
% xlabel('$t$');
% ylabel('$\tau_V$');
% grid on;

%%
% figure;
% plot(x/24,HI(:,1)); % HI(t,x)
% hold on;
% plot(x/24,HI(:,floor(length(x)/3)));
% plot(x/24,HI(:,floor(2*length(x)/3)));
% plot(x/24,HI(:,floor(length(x))));
% legend('$H_I(t,0)$','$H_I(t,1)$','$H_I(t,2)$','$H_I(t,3)$');
% xlabel('$t$');
% 
% figure;
% plot(x/24,VI(:,1)); % HI(t,x)
% hold on;
% plot(x/24,VI(:,floor(length(x)/3)));
% plot(x/24,VI(:,floor(2*length(x)/3)));
% plot(x/24,VI(:,floor(length(x))));
% legend('$V_I(t,0)$','$V_I(t,1)$','$V_I(t,2)$','$V_I(t,3)$');
% xlabel('$\tau_G$');

%%
toc