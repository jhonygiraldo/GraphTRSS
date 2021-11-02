clear all, close all, close all;
%%
line_width = 1.5;
marker_size = 6;
font_size = 21;
width = 680;
heigth = 460;
%%
load('condition_number_sob.mat');
%figure()
loglog(epsilon,condition_number_sob(:,1),'bs-','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
loglog(epsilon,condition_number_sob(:,2),'k*-','LineWidth',line_width,'MarkerSize',marker_size);
loglog(epsilon,condition_number_sob(:,3),'g-.','LineWidth',line_width,'MarkerSize',marker_size);
size_epsilon = length(epsilon);
semilogx(epsilon,condition_number_lap*ones(size_epsilon,1),'r-','LineWidth',line_width,'MarkerSize',marker_size);
lgd = legend({'$\kappa(\nabla_{\mathbf{z}}^2 f_{S}(\mathbf{z})),\beta=0.5$',...
    '$\kappa(\nabla_{\mathbf{z}}^2 f_{S}(\mathbf{z})),\beta=1$',...
    '$\kappa(\nabla_{\mathbf{z}}^2 f_{S}(\mathbf{z})),\beta=1.5$',...
    '$\kappa(\nabla_{\mathbf{z}}^2 f_{L}(\mathbf{z}))$'},'Location','northeast');
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
ylabel('$\kappa$','Interpreter','Latex');
xlabel('$\epsilon$','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex',...
    'XScale','log');
title('d) Condition Number','Interpreter','Latex');
xticks([1e-1, 1e0, 1e1]);
set(gcf,'Position',[100,100,width,heigth]);
axis tight;
ylim([250, 14000]);
grid on;
saveas(gcf,['condition_synthetic.svg']);