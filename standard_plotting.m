%% plot within-host dynamics
lineStyle = '-';
figure(1);
% subplot(1,3,1), plot(x/24,B,'-','LineWidth',3); hold on;
% title('$B(x)$','Interpreter','latex');
% hold on;
% subplot(1,3,2), plot(x/24,M,'-','LineWidth',3); hold on;
% title('$M(x)$','Interpreter','latex');
% xlabel('Age of infection (x) [days]');
%subplot(1,3,3), 
hold on;
plot(x/24,G,lineStyle,'LineWidth',3); 
xlim([0 600]);
ylim([1 (10.5)^5]);
title('Gametocyte abundance','Interpreter','latex');
xlabel('Time since infection (days)','Interpreter','latex');
yscale log;
length_infection_out = find(G>G_threshold,1,'last');
disp(['Length of infection within-host: ',num2str(h*length_infection_out/24),' days']);

figure(2);
hold on;
infection_level = h*sum(I,2);
plot(x/24,infection_level,lineStyle,'LineWidth',3); 
xlabel('Time since infection (days)','Interpreter','latex');
immune_activation = h*find(infection_level>P.IT,1,'first')/24;
if isempty(immune_activation)
    immune_activation = X_max;
end
xlim([0 600]);
ylim([1 (10.5)^6]);
yscale log;
title('iRBC abundance','Interpreter','latex');
%title('$\int I(x,\tau) \, d\tau$','Interpreter','latex');

% figure(3);
% hold on;
% plot(x/24,M,'-','LineWidth',3); 
% xlabel('Time since infection (days)','Interpreter','latex');
% xlim([0 400]);
% title('$M(x)$','Interpreter','latex');
% 
% figure(4);
% hold on;
% plot(x/24,B,'-','LineWidth',3); 
% xlabel('Time since infection (days)','Interpreter','latex');
% xlim([0 400]);
% title('$B(x)$','Interpreter','latex');

%legend('$B(x)$ uninfected red blood cells','$M(x)$ merozoites','$G(x)$ gametocytes');
%% Plot gametocytes and merozoite recruitment
lineStyle2 = ':';
figure(5);
hold on;
plot(x/24,G,lineStyle2,'Color',[0 0.4470 0.7410],'LineWidth',3);
plot(x/24,h*P.beta*sum(gamma_fun(tau,h).*I'),lineStyle2,'Color',[0.8500 0.3250 0.0980],'LineWidth',3);
%yline(10^2.4265,'--k','LineWidth',3); % sigmoid half maximum
%yline(10^6,'--k','LineWidth',3);
ylim([0.01 1000000]);
xline(immune_activation,lineStyle2,'Color','k','LineWidth',3);
yscale log;
xlim([0 100]);
legend('Gametocyte abundance','Merozoite recruitment','Interpreter','latex','Location','southeast');
xlabel('Time since infection (days)','Interpreter','latex');
%%
figure(6);
hold on;
yyaxis left
plot(x/24,P.sigma*(1-exp(-P.theta*A))*h.*sum(I,2),lineStyle2,'LineWidth',3); 
%xline(immune_activation,':k','LineWidth',3);
xlabel('Time since infection (days)','Interpreter','latex');
xlim([0 600]);
ylabel('$iRBC$ removal rate','Interpreter','latex');

yyaxis right
plot(x/24,(1-exp(-P.theta*A)),lineStyle2,'LineWidth',3); 
ylim([0 1]);
ylabel('Immune activation','Interpreter','latex');
%%
% figure(7);
% hold on;
% plot(x/24,I(:,1),'-','LineWidth',3); % I(x,tau)

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
% legend('$c = 1%$','$c = 4.4%$','$c = 10%$','Interpreter','latex','FontSize',35);

% Infectiousness plotting
figure(3);
plot(x/24,betaHV(G),lineStyle,'LineWidth',3); % Beta_HV(G(x))
hold on;
title('Host infectiousness','Interpreter','latex');
yline(betaHV(1000000000000000),'--k','LineWidth',3);
yline(betaHV(1000000000000000)/2,':k','LineWidth',3);
xlim([0 600]);
ylim([0 1]);
xlabel('Time since infection (days)','Interpreter','latex');
ylim([0 1]);

%% Plot cumulative infectiousness for given strategy
cum_inf1_time = h*cumtrapz(betaHV(G),1)/24;
figure(4);
hold on;
plot(x/24,cum_inf1_time,lineStyle,'Color',[0.9290    0.6940    0.1250],'LineWidth',3);
plot(x/24,h*cumtrapz(betaHV(10^10)*ones(length(x),1),1)/24,'--k','LineWidth',3); % upper bound strategy
xlim([0 600]);
ylim([0 300]);
title('Cumulative Infectiousness ($f_1$)','Interpreter','latex');
xlabel('Time since infection (days)','Interpreter','latex');
%legend('constant','$c = 25\%$','$c = 45\%$','$c = 60\%$','Interpreter','latex','FontSize',25);
% legend('$c = 4.4\%$','$c = 25\%$','$c = 45\%$','$c = 60\%$','Interpreter','latex','FontSize',25);
disp(['Strategy cumulative infectiousness within-host: ',num2str(cum_inf1_time(end))]);