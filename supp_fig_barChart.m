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
x = 0:1:4;
y = [267.86 274.08 282.31 283.73 284.25];
bar(x,y);
%text(0:length(y)-1,y,num2str(y'),'vert','bottom','horiz','center','FontSize',20); 
box off
axis tight;
xlabel('spline degree');
ylabel('cumulative infectiousness');
ylim([0 350]);
LimitsX = xlim; LimitsY = ylim;
title('C. fitness','FontWeight','Normal',...
    'HorizontalAlignment','left','position', [LimitsX(1), LimitsY(2)]);
set(gca,'TickDir','out');
box off;