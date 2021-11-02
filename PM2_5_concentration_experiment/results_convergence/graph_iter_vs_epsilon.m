clear all, close all, clc;
%%
load('iterations_vs_epsilon.mat');
load('epsilon_set.mat');
%%
line_width = 1.5;
marker_size = 6;
font_size = 21;
width = 680;
heigth = 460;
path_figures = 'figures_iter_vs_epsilon_vs_02/';
mkdir(path_figures);
%%
mean_errors = zeros(size(error_matrix_sobolev));
mean_iterations = zeros(size(repetitions_sobolev));
for i=1:size(error_matrix_sobolev,1)
    for j=1:size(error_matrix_sobolev,2)
        mean_errors(i,j) = mean(error_matrix_sobolev{i,j});
        mean_iterations(i,j) = mean(repetitions_sobolev{i,j});
    end
end
%%
yyaxis left
loglog(epsilon_set,mean(sqrt(mean_errors)),'bo--','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('RMSE','Interpreter','Latex');
yyaxis right
loglog(epsilon_set,mean(mean_iterations),'rs-','LineWidth',line_width,'MarkerSize',marker_size);
%semilogy(m,mean(repetitions_sobolev),'bs-','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('Iterations','Interpreter','Latex');
xlabel('$\epsilon$','Interpreter','Latex');
xlim([epsilon_set(1) epsilon_set(end)]);
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
xticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1]);
lgd = legend({'RMSE','Iterations'},'Location','northwest');
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
title('d) Variation $\epsilon$','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
grid on;
%saveas(gcf,[path_figures 'iterations_vs_epsilon.png']);
saveas(gcf,[path_figures 'iterations_vs_epsilon.svg']);