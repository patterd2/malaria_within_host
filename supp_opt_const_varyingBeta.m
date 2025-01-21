set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 25); % 25 for 1X3, 20 for 1X2
%set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%%
figure;
x = 12:0.5:17.5;
% beta = 12, 0.050390625 optimum
% beta = 13, 0.0646875   optimum
% beta = 14, 0.05375     optimum
% beta = 15, 0.050859375 optimum
% beta = 16, 0.04453125  optimum
% beta = 17, 0.03375     optimum
y = [0.050390625 0.059375 0.064686676025391 0.0559375...
     0.05375 0.054321350097656 0.050859375 0.047858459472656...
     0.04453125 0.039910949707031 0.03375 0.027104553222656];
scatter(x,100*y,200,'filled');
box off
axis tight;
xlabel('burst size');
ylabel('transmission investment');
ytickformat('percentage')
ylim([0 8]);
LimitsX = xlim; LimitsY = ylim;
title('A. optimal strategies','FontWeight','Normal',...
    'HorizontalAlignment','left','position', [LimitsX(1), LimitsY(2)]);
set(gca,'TickDir','out');
box off;