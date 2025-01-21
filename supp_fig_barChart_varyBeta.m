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
x = 1:3;
y = [121.52 187.30 283.73];
bar(x,y);
%text(0:length(y)-1,y,num2str(y'),'vert','bottom','horiz','center','FontSize',20); 
box off
axis tight;
xlabel('burst size');
ylim([0 300]);
yticks([0 100 200 300])
LimitsX = xlim; LimitsY = ylim;
title('cumulative infectiousness','FontWeight','Normal');
set(gca,'TickDir','out');
box off;