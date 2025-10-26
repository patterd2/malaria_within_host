%% Main file to run the full model
% (both within-host and human-vector systems)

tic
global P
set(0,'defaultTextFontName', 'Arial');
set(0,'defaultaxesfontsize', 20); % 25 for 1X3, 20 for 1X2 figures
%set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% Choose constant or nonconstant investment strategy

RUN_constant = 1;
RUN_nonconstant = 0;
RUN_degree = 3; % degree of the polynomial spline (1,2,3,4 - 3 baseline)

for i = 1:2
    %% numerical configuration
    X_max = 1000*24; % max time in days
    tau_max = 20*24; % max 20 days
    T_max = 200*24;
    xV_max = 20*24;
    G_threshold = 1; % gametocyte threshold to end infection (doesn't impact dynamics)
    h = 0.0625*2; % time/age step size in hours, same across all timescales

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
    if i == 1
        P.c = 0.044;
    else
        RUN_nonconstant = 1;
    end
    %% solve the within-host model
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

    % NB: ordering of independent variables is I(x,tau), IG(x,tau)


    %% Set the parasite investment strategy and run the within-host model

    if i == 1
        % vector CC stores the strategy (proportion investmented in onward transmission)
        CC = P.c*ones(1,nx); % set the baseline constant investment strategy
    else
        % Generate the strategy via a combination of the splines
        % * Note that grid must match the splines grid *
        if RUN_degree == 1
            temp1 = importdata('basisMatrixNoKnots_degree1_1000_0.125.txt'); % choose from spline files
            CC1 = temp1.data(:,1);
            CC2 = temp1.data(:,2);
            w1 = 0.094513570615327;
            w2 = -0.052554393112016;

            CC = min(1,max(0,w1*CC1 + w2*CC2))';
        elseif RUN_degree == 2
            temp1 = importdata('basisMatrixNoKnots_degree2_1000_0.125.txt'); % choose from spline files
            CC1 = temp1.data(:,1);
            CC2 = temp1.data(:,2);
            CC3 = temp1.data(:,3);
            w1 = 0.299711689528846;
            w2 = -0.374769344915881;
            w3 = 0.590730397030689;

            CC = min(1,max(0,w1*CC1 + w2*CC2 + w3*CC3))';
        elseif RUN_degree == 3
            temp1 = importdata('basisMatrixKnot75_degree3_1000_0.125.txt'); %

            CC1 = temp1.data(:,1);
            CC2 = temp1.data(:,2);
            CC3 = temp1.data(:,3);
            CC4 = temp1.data(:,4);

            % optimal weights with a knot at day 75, beta = 16
            w1 = -0.222178862091422;
            w2 = 0.512486123675708;
            w3 = -0.513683492605781;
            w4 = 0.478970226894363;

            CC = min(1,max(0,w1*CC1 + w2*CC2 + w3*CC3 + w4*CC4))';
        elseif RUN_degree == 4
            %temp1 = importdata('basisMatrixKnot75_degree4_1000_0.125.txt');
            temp1 = importdata('basisMatrixNoKnots_degree4_1000_0.125.txt'); % choose from spline files
            CC1 = temp1.data(:,1);
            CC2 = temp1.data(:,2);
            CC3 = temp1.data(:,3);
            CC4 = temp1.data(:,4);
            CC5 = temp1.data(:,5);
            w1 = 0.220470836603252; % weights with no knots
            w2 = 0.296449499531031;
            w3 = -0.790096277376987;
            w4 = 1.050393368686625;
            w5 = 0.005878963403393;
            % w1 = -0.221087421864143; % weights with knot at day 75
            % w2 = 0.512198995086951;
            % w3 = -0.510947054140931;
            % w4 = 0.470198208196176;
            % w5 = 0.027510246532859;

            CC = min(1,max(0,w1*CC1 + w2*CC2 + w3*CC3 + w4*CC4 + w5*CC5))';
        else
            CC = P.c*ones(1,nx); % run constant strat if no valid degree
        end
    end

    % uncomment the following lines of code to perturb the given strategy
    % sens_day = 580; % start time (in days) of the modification to strategy
    % sens_length = 20; % number of days for which strategy is modified
    % start_index = max(1,floor(sens_day*24/h));
    % CC(1,start_index:(start_index+floor(sens_length*24/h)) )...
    %     = 0.5*CC(1, start_index:(start_index+floor(sens_length*24/h)) );

    [B, M, I, IG, G, A] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

    %% solve the between-host/vector model
    % HS0 = 100; % scalar
    % HI0 = h*ones(1,nx); % HI(0,x), vector
    % VS0 = 100; % scalar
    % VI0 = zeros(1,nxV); % VI(t,xV) @ t = 0 (vector)
    %
    % [HS, HI, VS, VI] = human_vector_model(h, 0, X_max, T_max, xV_max, HS0, HI0, VS0, VI0, G);

    %% Plotting and sim info

    length_infection_out = find(G>G_threshold,1,'last');

    figure(1);
    hold on;
    if i == 1
        plot(x(1:length_infection_out)/24, 100*CC(1:length_infection_out),...
            'Color',[173/256 173/256 173/256],'LineWidth',3);
    else
        plot(x(1:length_infection_out)/24, 100*CC(1:length_infection_out),...
            'Color',[30/256 144/256 255/256],'LineWidth',3);
    end
    %scatter(x(length_infection_out)/24, 100*CC(length_infection_out),100,...
    %    [173/256 173/256 173/256],'diamond','filled')
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

    figure(2);
    hold on
    cum_inf1_time = h*cumtrapz(betaHV(G),1)/24;
    if i == 1
        plot(x(1:length_infection_out)/24,cum_inf1_time(1:length_infection_out),...
            '-','Color',[173/256 173/256 173/256],'LineWidth',3);
        temp_int = h*cumtrapz(betaHV(10^10)*ones(length(x),1),1)/24;
        plot(x/24,temp_int,'--k','LineWidth',3); % upper bound strategy
        scatter(x(length_infection_out)/24, cum_inf1_time(length_infection_out),...
            100,[173/256 173/256 173/256],'diamond','filled');
    else
        plot(x(1:length_infection_out)/24,cum_inf1_time(1:length_infection_out),...
            '-','Color',[30/256 144/256 255/256],'LineWidth',3);
        temp_int = h*cumtrapz(betaHV(10^10)*ones(length(x),1),1)/24;
        plot(x/24,temp_int,'--k','LineWidth',3); % upper bound strategy
        scatter(x(length_infection_out)/24, cum_inf1_time(length_infection_out),...
            100,[30/256 144/256 255/256],'diamond','filled');
    end
    xlim([0 600]);
    ylim([0 300]);
    xlabel('infection age x (days)');
    ylabel('cumulative infectiousness');
    LimitsX = xlim; LimitsY = ylim;
    title('B. fitness','FontWeight','Normal',...
        'HorizontalAlignment','left','position', [LimitsX(1), LimitsY(2)]);
    set(gca,'TickDir','out');
end

%%
figure(1);
hold on;
legend('constant strategy','age-varying strategy (cubic w/ knot day 75)');

toc