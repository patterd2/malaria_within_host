%% Main file to run the full model (both within-host and human-vector systems)
clc
format long
global P
tic

set(0,'defaultaxesfontsize', 20);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none')
set(0,'defaulttextinterpreter','latex');

%% numerical configuration
T_max = 8*24; % max time in days
P.T = T_max;
h = 0.05; % time/age since infection, etc. step size in hours;
x = (0:h:T_max)';
nx = length(x);

% set model parameters via the baseline file (contains global variables)
baseline_parameter_set;

%% solve the within-host model
% initially there are no merozoites or (developing/mature) gametocytes
B0 = P.Bstar; % scalar, nonzero
M0 = 0; % scalar, zero
I0 = h*ones(1,nx); % I(0,tau), should be nonzero, preserve integral (total number of infected)
IG0 = zeros(1,nx); % IG(0,tau)
G0 = 0; % scalar, zero
A0 = 0; % scalar, zero

% NB: ordering of independent variables is I(x,tau), IG(x,tau)

[B, M, I, IG, G, A] = within_host_model(h, 0, T_max, B0, M0, I0, IG0, G0, A0);
%% solve the between-host/vector model
HS0 = 100; % scalar
HI0 = ones(1,nx); % vector
VS0 = 100; % scalar
VI0 = zeros(1,nx); % vector

[HS, HI, VS, VI] = human_vector_model(h, 0, T_max, HS0, HI0, VS0, VI0, G);
%% plot within-host dynamics
figure;
plot(x/24,B,'LineWidth',3);
hold on;
plot(x/24,M,'LineWidth',3);
plot(x/24,G,'LineWidth',3);
%plot(x/24,A,'LineWidth',2);
xlabel('Time since infection (days)');
legend('$B(x)$ uninfected red blood cells','$M(x)$ merozoites','$G(x)$ gametocytes');
grid on;

%%
figure;
plot(x/24,I(:,1)); % I(x,tau)
hold on;
plot(x/24,I(:,floor(length(x)/6)));
plot(x/24,I(:,floor(length(x)/3)));
plot(x/24,I(:,floor(2*length(x)/3)));
plot(x/24,I(:,floor(length(x))));
title('Infection dynamics (asexual stage): $I(x,\tau)$');
xlabel('Time since infection (days)');
legend('$I(x,0)$','$I(x,0.5)$','$I(x,1)$','$I(x,2)$','$I(x,3)$');
axis tight;
grid on;

%% 
figure;
plot(x/24,sum(I,2)); % I(x,tau)
title('Infection dynamics (asexual stage): $\int I(x,\tau) \, d\tau$');
xlabel('Time since infection (days)');
axis tight;
grid on;

%%

figure;
plot(x/24,IG(:,1)); % I(x,tau)
hold on;
plot(x/24,IG(:,floor(length(x)/6)));
plot(x/24,IG(:,floor(length(x)/3)));
plot(x/24,IG(:,floor(2*length(x)/3)));
plot(x/24,IG(:,floor(length(x))));
title('Infection dynamics (sexual stage): $I_G(x,\tau)$');
xlabel('Time since infection (days)');
legend('$I_G(x,0)$','$I_G(x,0.5)$','$I_G(x,1)$','$I_G(x,2)$','$I_G(x,3)$');
axis tight;
grid on;

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

% figure;
% imagesc(x/24,x/24,HI);
% title('$H_I(t,x)$');
% set(gca,'YDir','normal');
% colorbar;
% clim([min(min(HI)) max(max(HI))+1]);
% xlabel('$t$');
% ylabel('$x$');
% grid on;

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
toc;