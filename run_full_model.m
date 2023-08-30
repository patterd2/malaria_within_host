%% Main file to run the full model (both within-host and human-vector systems)
clc
format long
global P
tic
%% numerical configuration
T_max = 3*24; % max time in days
P.T = T_max;
h = 0.01; % time/age since infection, etc. step size in hours;
x = (0:h:T_max)';
nx = length(x);

% set model parameters via the baseline file (contains global variables)
baseline_parameter_set;

%% solve the within-host model
% initially there are no merozoites or (developing/mature) gametocytes
B0 = 100; % scalar, nonzero
M0 = 0; % scalar, zero
I0 = ones(1,nx); % vector, should be nonzero
IG0 = zeros(1,nx); % vector, zero
G0 = 0; % scalar, zero
A0 = 0; % scalar, zero

[B, M, I, IG, G, A] = within_host_model(h, 0, T_max, B0, M0, I0, IG0, G0, A0);
%% solve the between-host/vector model
HS0 = 100; % scalar
HI0 = ones(1,nx); % vector
VS0 = 100; % scalar
VI0 = zeros(1,nx); % vector

[HS, HI, VS, VI] = human_vector_model(h, 0, T_max, HS0, HI0, VS0, VI0, G);
%% plot within-host dynamics
figure_setups;
plot(x/24,B,'LineWidth',2);
hold on;
plot(x/24,M,'LineWidth',2);
plot(x/24,G,'LineWidth',2);
%plot(x/24,A,'LineWidth',2);
xlabel('Time since infection (days)');
legend('(uninfected) red blood cells','merozoites','gametocytes');
grid on;

figure;
imagesc(x/24,x/24,I);
title('$I(x,\tau)$');
set(gca,'YDir','normal');
colorbar;
clim([min(min(I)) max(max(I))+1]);
xlabel('$x$');
ylabel('$\tau$');
grid on;

figure;
imagesc(x/24,x/24,IG);
title('$I_G(x,\tau_G)$');
set(gca,'YDir','normal');
colorbar;
clim([min(min(IG)) max(max(IG))+1]);
xlabel('$x$');
ylabel('$\tau_G$');
grid on;

%% plot human-vector dynamics
figure_setups;
plot(x/24,HS,'LineWidth',2);
hold on;
plot(x/24,VS,'LineWidth',2);
xlabel('Time (days)');
legend('susceptible humans','susceptible vectors');
grid on;

figure;
imagesc(x/24,x/24,HI);
title('$H_I(t,x)$');
set(gca,'YDir','normal');
colorbar;
clim([min(min(HI)) max(max(HI))+1]);
xlabel('$t$');
ylabel('$x$');
grid on;

figure;
imagesc(x/24,x/24,VI);
title('$V_I(t,\tau_V)$');
set(gca,'YDir','normal');
colorbar;
clim([min(min(VI)) max(max(VI))+1]);
xlabel('$t$');
ylabel('$\tau_V$');
grid on;

%%
toc;