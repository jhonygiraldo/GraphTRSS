clear all, close all, clc;
%%
load('../graph_construction/full_graph.mat');
load('../covid_19_new_cases.mat');
G = gsp_compute_fourier_basis(G);
[N,T] = size(Data);
Dh = temporal_difference_operator(T);
XDh = Data*Dh;
X_hat = G.U'*XDh;
%%
cutoff_energy = 0.9;
X_hat = X_hat.^2;
total_energy = sum(X_hat);
cutoff_energy = total_energy*cutoff_energy;
X_hat_cum = cumsum(X_hat);
bandwidths = zeros(T-1,1);
for i=1:T-1
    bandwidths(i) = max(find(X_hat_cum(:,i) <= cutoff_energy(i)));
end
mean(bandwidths)/N
%%
line_width = 1.5;
marker_size = 6;
font_size = 21;
width = 680;
heigth = 460;
path_figures = 'figures/';
mkdir(path_figures);
%%
indx_to_show = [173,233,293,104];
for i=1:length(indx_to_show)
    figure()
    plot(X_hat(:,indx_to_show(i)),'LineWidth',line_width,'MarkerSize',marker_size);
    ylabel(['$(\hat{\mathbf{x}}_{',num2str(indx_to_show(i)+1),...
        '}-\hat{\mathbf{x}}_{',num2str(indx_to_show(i)),'})^{(2)}$'],'Interpreter','Latex');
    xlabel('Eigenvalue Index','Interpreter','Latex');
    title(['Bandwidth = ',num2str(bandwidths(indx_to_show(i)))],'Interpreter','Latex');
    get(gca);
    set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
    set(gcf,'Position',[100,100,width,heigth]);
    grid on;
    xlim([0,N]);
    saveas(gcf,[path_figures 'bandwidth' num2str(i) '.svg']);
end