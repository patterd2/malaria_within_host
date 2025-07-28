%% Script to perform strategy optimization for nonconstant parasite investment
tic

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 25 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

% NB Make sure parameters agree with withinhost_model_optimization.m
baseline_parameter_set;
X_max = 1000*24; % max time in days
h = 0.0625*2; % time/age step size in hours, same across all timescales
x = (0:h:X_max)';
N = 1; % number of initial conditions to try
max_cum_inf = zeros(N,1);
save_strats = zeros(N,4); % matrix to save optimal weights
save_init = zeros(N,4); %

tic;
for i = 1:N
    % the dimension of v will set the strategy space
    % e.g. length(v) = 3 selects for quadratic splines, etc.
    %v = [-0.222178862091422+0.5*(rand()-0.5) 0.512486123675708+0.5*(rand()-0.5)...
    %     -0.513683492605781+0.5*(rand()-0.5) 0.478970226894363+0.5*(rand()-0.5)];
    %v = [rand()-0.5 rand()-0.5 rand()-0.5 rand()-0.5];
    v = [-0.473309466837109 0.731820507032210 -1.358692478967953 5.734801386903090];
    options = optimset('Display','iter','MaxIter',100);
    %options = optimset('Display','iter','MaxIter',50,'PlotFcns','optimplotfval','TolX',1e-7);
    %options = optimset('MaxIter',100);
    [a, funmax] = fminsearch(@withinhost_model_optimization,v,options);
    %toc
    max_cum_inf(i) = -funmax;
    save_strats(i,1:length(a)) = a;

    output_string = sprintf('Run %d: Optimal cumulative infectiousness %f', i, -funmax);
    disp(output_string); % output optimal fitness value
end
toc; % end timing when optimization is complete
%% Plot all strategies
for i = 1:N
    a = save_strats(i,:);
    if isscalar(a)
        CC = a*ones(1,length(x));
    elseif length(a) == 2
        temp1 = importdata('basisMatrixNoKnots_degree1_1000_0.125.txt'); % choose linear splines
        %temp1 = importdata('basisMatrixKnot.txt'); % choose from spline files
        CC1 = temp1.data(:,1);
        CC2 = temp1.data(:,2);
        CC = min(1,max(0,a(1)*CC1 + a(2)*CC2));
        CC_init = min(1,max(0,v(1)*CC1 + v(2)*CC2));
    elseif length(a) == 3
        temp1 = importdata('basisMatrixNoKnots_degree2_1000_0.125.txt'); % choose quadratic splines
        CC1 = temp1.data(:,1);
        CC2 = temp1.data(:,2);
        CC3 = temp1.data(:,3);
        CC = min(1,max(0,a(1)*CC1 + a(2)*CC2 + a(3)*CC3));
        CC_init = min(1,max(0,v(1)*CC1 + v(2)*CC2 + v(3)*CC3));
    elseif length(a) == 4
        %temp1 = importdata('basisMatrixKnot125_degree3_1000_0.125.txt');
        %temp1 = importdata('basisMatrixKnot25_degree3_1000_0.125.txt'); % with knot
        temp1 = importdata('basisMatrixKnot75_degree3_1000_0.125.txt'); % with knot
        %temp1 = importdata('basisMatrixNoKnots_1000_0.125.txt'); % no knot
        CC1 = temp1.data(:,1);
        CC2 = temp1.data(:,2);
        CC3 = temp1.data(:,3);
        CC4 = temp1.data(:,4);
        CC = min(1,max(0,a(1)*CC1 + a(2)*CC2 +...
            a(3)*CC3 + a(4)*CC4));% + a(5)*CC5);
        CC_init = min(1,max(0,v(1)*CC1 + v(2)*CC2 +...
            v(3)*CC3 + v(4)*CC4)); % + v(5)*CC5);
    elseif length(a) == 5
        temp1 = importdata('basisMatrixKnot75_degree4_1000_0.125.txt'); % with knot
        %temp1 = importdata('basisMatrixNoKnots_degree4_1000_0.125.txt'); % no knot
        CC1 = temp1.data(:,1);
        CC2 = temp1.data(:,2);
        CC3 = temp1.data(:,3);
        CC4 = temp1.data(:,4);
        CC5 = temp1.data(:,5);
        CC = min(1,max(0,a(1)*CC1 + a(2)*CC2 +...
            a(3)*CC3 + a(4)*CC4));% + a(5)*CC5);
        CC_init = min(1,max(0,v(1)*CC1 + v(2)*CC2 +...
            v(3)*CC3 + v(4)*CC4+ v(5)*CC5));
    end
    % Plot optimal strategy
    figure(1);
    hold on;
    plot(x/24,100*CC,'LineWidth',3);
    ylim([0 100]);
    xlim([0 600]);
    xlabel('infection age x (days)');
    ytickformat('percentage');
end

%% plot the optimal strategy versus time

[~, opt_strat] = max(max_cum_inf); % identify the optimum
a = save_strats(opt_strat,:);

if isscalar(a)
    CC = a*ones(1,length(x));
elseif length(a) == 2
    temp1 = importdata('basisMatrixNoKnots_degree1_1000_0.125.txt'); % choose linear splines
    %temp1 = importdata('basisMatrixKnot.txt'); % choose from spline files
    CC1 = temp1.data(:,1);
    CC2 = temp1.data(:,2);
    a = save_strats(opt_strat, :);
    CC = min(1,max(0,a(1)*CC1 + a(2)*CC2));
    CC_init = min(1,max(0,v(1)*CC1 + v(2)*CC2));
elseif length(a) == 3
    temp1 = importdata('basisMatrixNoKnots_degree2_1000_0.125.txt'); % choose quadratic splines
    CC1 = temp1.data(:,1);
    CC2 = temp1.data(:,2);
    CC3 = temp1.data(:,3);
    a = save_strats(opt_strat, :);
    CC = min(1,max(0,a(1)*CC1 + a(2)*CC2 + a(3)*CC3));
    CC_init = min(1,max(0,v(1)*CC1 + v(2)*CC2 + v(3)*CC3));
elseif length(a) == 4
    %temp1 = importdata('basisMatrixKnot125_degree3_1000_0.125.txt');
    %temp1 = importdata('basisMatrixKnot25_degree3_1000_0.125.txt'); % with knot
    temp1 = importdata('basisMatrixKnot75_degree3_1000_0.125.txt'); % with knot
    %temp1 = importdata('basisMatrixNoKnots_1000_0.125.txt'); % no knot
    CC1 = temp1.data(:,1);
    CC2 = temp1.data(:,2);
    CC3 = temp1.data(:,3);
    CC4 = temp1.data(:,4);
    a = save_strats(opt_strat, :);
    CC = min(1,max(0,a(1)*CC1 + a(2)*CC2 +...
        a(3)*CC3 + a(4)*CC4));% + a(5)*CC5);
    CC_init = min(1,max(0,v(1)*CC1 + v(2)*CC2 +...
        v(3)*CC3 + v(4)*CC4)); % + v(5)*CC5);
elseif length(a) == 5
    temp1 = importdata('basisMatrixKnot75_degree4_1000_0.125.txt'); % with knot
    %temp1 = importdata('basisMatrixNoKnots_degree4_1000_0.125.txt'); % no knot
    CC1 = temp1.data(:,1);
    CC2 = temp1.data(:,2);
    CC3 = temp1.data(:,3);
    CC4 = temp1.data(:,4);
    CC5 = temp1.data(:,5);
    a = save_strats(opt_strat, :);
    CC = min(1,max(0,a(1)*CC1 + a(2)*CC2 +...
        a(3)*CC3 + a(4)*CC4));% + a(5)*CC5);
    CC_init = min(1,max(0,v(1)*CC1 + v(2)*CC2 +...
        v(3)*CC3 + v(4)*CC4+ v(5)*CC5));
end
% Plot optimal strategy
figure(1);
hold on;
plot(x/24,100*CC,'--k','LineWidth',3);
ylim([0 100]);
xlim([0 600]);
xlabel('infection age x (days)');
ytickformat('percentage');