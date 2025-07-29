a = [-0.430376295532938 0.705412349737878 -0.897850168548958 3.822906252730528];

% uncomment to choose a spline basis
%temp1 = importdata('basisMatrixKnot125_degree3_1000_0.125.txt');
%temp1 = importdata('basisMatrixKnot25_degree3_1000_0.125.txt'); % with knot
temp1 = importdata('basisMatrixKnot75_degree3_1000_0.125.txt'); % with knot
%temp1 = importdata('basisMatrixNoKnots_1000_0.125.txt'); % no knot
CC1 = temp1.data(:,1);
CC2 = temp1.data(:,2);
CC3 = temp1.data(:,3);
CC4 = temp1.data(:,4);
CC = min(1,max(0,a(1)*CC1 + a(2)*CC2 +...
    a(3)*CC3 + a(4)*CC4));
disp(withinhost_model_optimization(a)); % print the fitness of the strategy
figure(1);
hold on;
plot(x/24,100*CC,'LineWidth',3);
ylim([0 100]);
xlim([0 600]);
xlabel('infection age x (days)');
ytickformat('percentage');