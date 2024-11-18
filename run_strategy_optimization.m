%% Script to perform strategy optimization for nonconstant parasite investment

% NB Make sure parameters agree with withinhost_model_optimization.m
baseline_parameter_set;
X_max = 280*24; % max time in days
h = 0.25; % time/age step size in hours, same across all timescales
x = (0:h:X_max)';
N = 1; % number of initial conditions to try
max_cum_inf = zeros(N,1);
save_strats = zeros(N,4); % matrix to save optimal weights
save_init = zeros(N,4); %

tic;
for i = 1:N
    %tic
    %v = [-0.223 0.506 -0.282 0.153]; % initial weights for cubic splines
    v = [(rand()-0.5) (rand()-0.5) (rand()-0.5) (rand()-0.5)];
    options = optimset('Display','iter','MaxIter',200);
    %options = optimset('MaxIter',100);
    [a, funmax] = fminsearch(@withinhost_model_optimization,v,options);
    %toc
    max_cum_inf(i) = -funmax;
    save_strats(i,:) = a;

    output_string = sprintf('Run %d: Optimal cumulative infectiousness %f', i, -funmax);
    disp(output_string); % output optimal weights
end
toc;
%% plot the optimal strategy versus time
[~, opt_strat] = max(max_cum_inf); % identify the optimum 

% July 9th optimal plot on Overleaf:
% w1 = 0.199944393818933; w2 = 0.197075562766502; w3 = -0.453310918950100;
% w4 = 0.714969288272976;
% Strategy cumulative infectiousness within-host: 297.9466

% generate the strategy
temp1 = importdata('basisMatrixNoKnots_1000_0.5.txt'); % choose from spline files
%temp1 = importdata('basisMatrixKnot.txt'); % choose from spline files
CC1 = temp1.data(:,1);
CC2 = temp1.data(:,2);
CC3 = temp1.data(:,3);
CC4 = temp1.data(:,4);
%CC5 = temp1.data(:,5);
a = save_strats(opt_strat, :);
CC = min(1,max(0,a(1)*CC1 + a(2)*CC2 +...
    a(3)*CC3 + a(4)*CC4));% + a(5)*CC5);
CC_init = min(1,max(0,v(1)*CC1 + v(2)*CC2 +...
    v(3)*CC3 + v(4)*CC4)); % + v(5)*CC5);

figure;
hold on;
plot(x/24,100*P.c*ones(1,length(x)),'LineWidth',3);
plot(x/24,100*CC,'LineWidth',3);
%plot(x/24,100*CC_init,'--','LineWidth',3);
ylim([0 100]);
xlabel('Time since infection (days)','Interpreter','latex');
ytickformat('percentage');
legend('Optimal Strategy');

