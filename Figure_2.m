%% Code to produce Figure 2

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

for i = 1:2
    % set model parameters via the baseline file (contains global variables)
    baseline_parameter_set;
    P.c = 0.05;
    if i == 2
        P.sigma = 0;
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

    CC = P.c*ones(1,nx); % set strategy

    [B, M, I, IG, G, A] = within_host_model(h, 0, X_max, tau_max, B0, M0, I0, IG0, G0, A0, CC);

    %% Plotting and sim info
    
    if i == 1
        lineStyle = '-';
    else
        lineStyle = '--';
    end

    figure(3);
    hold on;
    infection_level = h*sum(I,2);
    plot(x/24,infection_level,lineStyle,'color',[0.4219 0.3867 0.6758],'LineWidth',3);
    xlabel('infection age x (days)');
    xlim([0 540]);
    xticks([0 120 240 360 480])
    ylim([1 (10.5)^6]);
    ylabel('abundance per \muL','Interpreter','tex');
    yscale log;
    %title('iRBC abundance','FontWeight','normal');
    set(gca,'TickDir','out');

    plot(x/24,G,lineStyle,'color',[205/256 101/256 146/256],'LineWidth',3);

    set(gca,'TickDir','out');
    title('C. parasite dynamics','FontWeight','normal');
    legend('iRBC','gametocyte','Location','south');
    legend('boxoff');

    figure(1);
    hold on;
    immune_removal = 100*(1-exp(-P.theta*A));
    if i == 1
        plot(x/24,immune_removal,lineStyle,'color',[250/256 196/256 131/256],'LineWidth',3);
    else
        plot(x/24,zeros(length(x),1),lineStyle,'color',[250/256 196/256 131/256],'LineWidth',3);
    end
    xlabel('infection age x (days)');
    xlim([0 540]);
    xticks([0 120 240 360 480])
    ylim([-1 100]);
    ylabel('% maximal immune removal');

    set(gca,'TickDir','out');
    title('B. immune dynamics','FontWeight','normal');
    if i == 2
        legend('immunity','no immunity','Location','south');
        legend('boxoff');
    end

    figure(2);
    hold on;
    % find(infection_level>10,1,'first')
    dA = diff(A)/h;
    if i == 1
        plot(infection_level(1:2039),dA(1:2039),lineStyle,'color',[250/256 196/256 131/256],'LineWidth',3);
    else
        plot(infection_level(1:2039),zeros(2039,1),lineStyle,'color',[250/256 196/256 131/256],'LineWidth',3);
    end
    xlabel('iRBCs per \muL','Interpreter','tex');
    xlim([0 10]);
    %xticks([0 120 240 360 480])
    ylim([-0.01 1]);
    ylabel('dA(x)/dx');

    set(gca,'TickDir','out');
    title('B. immune upregulation (inset)','FontWeight','normal');
    if i == 2
        legend('immunity','no immunity','Location','south');
        legend('boxoff');
    end


end
