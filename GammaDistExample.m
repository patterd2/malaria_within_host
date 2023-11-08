%% Creating and plotting gamma distribution

x = linspace(0,15,1000); % x-values to consider

sc = 10; % number of identical states (more states gives narrower distribution)
a1 = 47; % shape parameter
b1 = 1/4.47; % scale parameter
mean1 = a1*b1; % mean = shape * scale

% Create gamma distribution pdf
y1 = gampdf(x,a1,b1);
% Create gamma distribution cdf
y1cdf = gamcdf(x,a1,b1);

%% Redo with different shape parameter (to shift mean)

% sc2 = sc; % number of identical states
% a2 = 47; % shape parameter
% b2 = 4.47; % scale parameter
% mean2 = a2*b2; % mean = shape * scale
% 
% % Create gamma distribution pdf
% y2 = gampdf(x,a2,b2);
% % Create gamma distribution cdf
% y2cdf = gamcdf(x,a2,b2);

%% Plot gamma distribution pdf

figure(10)
plot(x,y1,'b','linewidth',2)
hold on
plot(x,y2,'r','linewidth',2)
plot([mean1 mean1],[0 max([y1 y2])*1.1],'b')
plot([mean2 mean2],[0 max([y1 y2])*1.1],'r')
hold off
set(gca,'fontsize',14)
title('Probability Distribution Function')

%% Plot gamma distribution cdf

figure(11)
plot(x,y1cdf,'b','linewidth',2)
hold on
plot(x,y2cdf,'r','linewidth',2)
plot([mean1 mean1],[0 1],'b')
plot([mean2 mean2],[0 1],'r')
hold off
grid on
set(gca,'fontsize',14)
title('Cumulative Distribution Function')

