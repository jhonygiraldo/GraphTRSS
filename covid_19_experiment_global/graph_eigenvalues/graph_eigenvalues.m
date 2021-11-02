clear all, close all, clc;
%%
line_width = 1.5;
marker_size = 6;
font_size = 21;
width = 680;
heigth = 410;
%%
load('../graph_construction/full_graph.mat');
[V,D] = eig(full(G.L));
eigenvalues = sort(diag(D));
beta = [0.1, 0.3, 0.6, 1, 2, 5];
LegendInfo{1} = ['$\beta=',num2str(beta(1)),'$']
plot(eigenvalues.^beta(1)/(max(eigenvalues.^beta(1))),'LineWidth',...
    line_width,'MarkerSize',marker_size);%,'DisplayName',txt);
hold on;
for i=2:length(beta)
    LegendInfo{i} = ['$\beta=',num2str(beta(i)),'$']
    plot(eigenvalues.^beta(i)/(max(eigenvalues.^beta(i))),'LineWidth',...
        line_width,'MarkerSize',marker_size);%,'DisplayName',txt);
end
ylabel('Magnitude of Eigenvalue $\lambda_i$','Interpreter','Latex');
xlabel('Index of Eigenvalue $i$','Interpreter','Latex');
lgd = legend(LegendInfo,'Location','northeast');
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
title(['Eigenvalues Penalization for COVID-19'],'Interpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
axis tight;
grid on;
saveas(gcf,['eigenvalue_penalization.svg']);