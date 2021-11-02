clear all, close all, close all;
%%
line_width = 1.5;
marker_size = 6;
font_size = 21;
width = 680;
heigth = 460;
%%
load('condition_number_sob.mat');
semilogx(epsilon,condition_number_sob(:,1),'bs-','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
size_epsilon = length(epsilon);
semilogx(epsilon,condition_number_lap*ones(size_epsilon,1),'r-','LineWidth',line_width,'MarkerSize',marker_size);
lgd = legend({'$\kappa(\nabla_{\mathbf{z}}^2 f_{S}(\mathbf{z}))$',...
    '$\kappa(\nabla_{\mathbf{z}}^2 f_{L}(\mathbf{z}))$'},'Location','northeast');
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
ylabel('Condition Number $\kappa$','Interpreter','Latex');
xlabel('$\epsilon$','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex',...
    'XScale','log');
title('a) COVID-19 Global, $\beta=1$','Interpreter','Latex');
xticks([1e-4, 1e-2, 1e0, 1e2, 1e4]);
set(gcf,'Position',[100,100,width,heigth]);
axis tight;
grid on;
saveas(gcf,['condition_covid_global.svg']);