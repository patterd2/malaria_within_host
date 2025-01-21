%% plot within-host dynamics
lineStyle = '-';
figure(1);
hold on;
infection_level = h*sum(I,2);
plot(x/24,infection_level,lineStyle,'LineWidth',3); 
xlabel('infection age x (days)');
immune_activation = h*find(infection_level>P.IT,1,'first')/24;
if isempty(immune_activation)
    immune_activation = X_max;
end
xlim([0 800]);
ylim([1 (10.5)^6]);
ylabel('abundance');
yscale log;
%title('iRBC abundance','FontWeight','normal');
set(gca,'TickDir','out');

plot(x/24,G,lineStyle,'LineWidth',3); 
length_infection_out = find(G>G_threshold,1,'last');
disp(['Length of infection within-host: ',num2str(h*length_infection_out/24),' days']);
set(gca,'TickDir','out');
title('parasite dynamics','FontWeight','normal');
legend('iRBC','gametocyte');
legend('boxoff');

%%
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
%lineStyle2 = '-';
% figure(5);
% hold on;
% plot(x/24,G,lineStyle2,'Color',[0 0.4470 0.7410],'LineWidth',3);
% plot(x/24,h*P.beta*sum(gamma_fun(tau,h).*I'),lineStyle2,'Color',[0.8500 0.3250 0.0980],'LineWidth',3);
% yline(10^2.4265,'--k','LineWidth',3); % sigmoid half maximum
% yline(10^6,'--k','LineWidth',3);
% ylim([0.01 1000000]);
% %xline(immune_activation,lineStyle2,'Color','k','LineWidth',3);
% yscale log;
% xlim([0 100]);
% legend('Gametocyte abundance','Merozoite recruitment','Interpreter','latex','Location','southeast');
% xlabel('Time since infection (days)','Interpreter','latex');
%%
% figure(6);
% hold on;
% yyaxis left
% plot(x/24,P.sigma*(1-exp(-P.theta*A))*h.*sum(I,2),'.-','LineWidth',3); 
% %xline(immune_activation,':k','LineWidth',3);
% xlabel('infection age x (days)');
% xlim([0 600]);
% ylabel('iRBC removal rate');
% set(gca,'TickDir','out');
% 
% yyaxis right
% plot(x/24,(1-exp(-P.theta*A)),'.-','LineWidth',3); 
% ylim([0 1]);
% ylabel('Immune activation');
% set(gca,'TickDir','out');
%% Plot A and dA/dx
% figure;
% plot(x(1:end-1)/24,diff(A)/h,'LineWidth',3);
% xlabel('days');
% ylabel('dA/dx');
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

%% Strategy plotting
figure(11);
hold on;
plot(x(1:length_infection_out)/24, 100*CC(1:length_infection_out),'LineWidth',3);
%plot(x(length_infection_out:end)/24, 100*CC(length_infection_out:end),'--','LineWidth',3);
scatter(x(length_infection_out)/24, 100*CC(length_infection_out),100,'k','diamond','filled')
axis tight;
xlabel('infection age x (days)');
ylabel('transmission investment');
ytickformat('percentage');
ylim([0 50]);
xlim([0 600]);
LimitsX = xlim; LimitsY = ylim;
title('A. optimal strategies','FontWeight','Normal',...
    'HorizontalAlignment','left','position', [LimitsX(1), LimitsY(2)]);
set(gca,'TickDir','out');
box off;

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
% figure(3);
% hold on;
% plot(x(1:length_infection_out)/24, betaHV(G(1:length_infection_out)),'LineWidth',3);
% %plot(x(length_infection_out:end)/24, betaHV(G(length_infection_out:end)),'--','LineWidth',3);
% scatter(x(length_infection_out)/24, betaHV(G(length_infection_out)),100,'k','diamond','filled')
% yline(betaHV(1000000000000000),'--k','LineWidth',3);
% %yline(betaHV(1000000000000000)/2,':k','LineWidth',3);
% xlim([0 1000]);
% ylim([0 1]);
% xlabel('infection age x (days)');
% ylabel('prop. mosquitoes infected');
% ylim([0 1]);
% %LimitsX = xlim; LimitsY = ylim;
% %title('B. host infectiousness','FontWeight','Normal',...
% %    'HorizontalAlignment','left','position', [LimitsX(1), LimitsY(2)]);
% set(gca,'TickDir','out');
% box off;

%% Merozoite reproduction plotting

% V1
% P1 = h*trapz(I(:,:).*repmat(gamma_fun(tau,h),1,length(x))',2);
% P2 = h*trapz(I(:,:),2)*P.mu;
% P3 = h*trapz(I(:,:),2).*P.sigma.*(1-exp(-P.theta*A(:)));
% 
% R_M = (1-P.c).*P.beta*(P.p*B(:,:)./(P.p*B(:,:)+P.muM))...
%     .*(P1./(P1 + P2 + P3)); % slow calculation

% P1 = trapz(I(:,:).*repmat(gamma_fun(tau,h),1,length(x))',2)./trapz(I(:,:),2);
% P2 = P.mu;
% P3 = P.sigma.*(1-exp(-P.theta*A(:)));
% 
% R_M = (1-P.c).*P.beta*(P.p*B(:,:)./(P.p*B(:,:)+P.muM))...
%     .*(P1./(P1 + P2 + P3)); % slow calculation
% 
% figure(13);
% hold on;
% plot(x(:)/24,R_M,'.-','LineWidth',3);
% scatter(x(length_infection_out)/24, R_M(length_infection_out),100,'k','diamond','filled');
% xlabel('age of infection (days)');
% title('Effective Merozoite Number','FontWeight','Normal');
% axis tight;

%% Plot cumulative infectiousness for given strategy
cum_inf1_time = h*cumtrapz(betaHV(G),1)/24;
figure(4);
hold on;
plot(x(1:length_infection_out)/24,cum_inf1_time(1:length_infection_out),'-','LineWidth',3);
temp_int = h*cumtrapz(betaHV(10^10)*ones(length(x),1),1)/24;
plot(x/24,temp_int,'--','LineWidth',3); % upper bound strategy
scatter(x(length_infection_out)/24, cum_inf1_time(length_infection_out),100,'k','diamond','filled')
xlim([0 600]);
ylim([0 300]);
xlabel('infection age x (days)');
ylabel('cumulative infectiousness');
LimitsX = xlim; LimitsY = ylim;
title('B. host infectiousness','FontWeight','Normal',...
    'HorizontalAlignment','left','position', [LimitsX(1), LimitsY(2)]);
set(gca,'TickDir','out');
disp(['Strategy cumulative infectiousness within-host: ',num2str(cum_inf1_time(end))]);

%% Plot the rate at which max infectiousness is achieved
cum_inf1_time = h*cumtrapz(betaHV(G),1)/24;
figure(14);
hold on;
plot(x(1:length_infection_out)/24,100*cum_inf1_time(1:length_infection_out)/cum_inf1_time(end),'-','LineWidth',3);
%temp_int = h*cumtrapz(betaHV(10^10)*ones(length(x),1),1)/24;
%plot(x/24,temp_int,'--','LineWidth',3); % upper bound strategy
scatter(x(length_infection_out)/24, 100*cum_inf1_time(length_infection_out)/cum_inf1_time(end),100,'k','diamond','filled')
ytickformat('percentage');
xlim([0 600]);
ylim([0 100]);
xlabel('infection age x (days)');
ylabel('prop. max infectiousness');
LimitsX = xlim; LimitsY = ylim;
title('B. host infectiousness','FontWeight','Normal',...
    'HorizontalAlignment','left','position', [LimitsX(1), LimitsY(2)]);
set(gca,'TickDir','out');
disp(['Strategy cumulative infectiousness within-host: ',num2str(cum_inf1_time(end))]);