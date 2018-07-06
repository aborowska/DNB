ff = figure(123);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.45 0.95]);
set(gcf,'defaulttextinterpreter','latex');

for ii = 1:7
    subplot(4,2,ii)
    MAT = squeeze(ChainsDNB_GT(:,ii,:));
    boxplot(MAT)
    hold on
    plot(get(gca,'xlim'),[param_true(ii),param_true(ii)],'k')
    hold off
    title(param{ii})
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Simulation no.','Interpreter','latex')
%     ll = legend('True value');
%     set(ll,'Interpreter','latex','Location','SouthEast')
end

set(gcf,'PaperPositionMode','auto');
plotname = 'DNB_sim_barplot.png';
print(ff,plotname,'-dpng','-r0')
plotname = 'DNB_sim_barplot.eps';
print(ff,plotname,'-depsc','-r0')
