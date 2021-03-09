clear all;
close all;
data = dlmread('data_error_grad.txt','',1,0);
n_dofs = data(:,1);
P2UMF = data(:,2);
P2CG_1em16 = data(:,3);
P2CG_1em10 = data(:,4);
P2CG_1em4 = data(:,5);

h = figure;
% loglog(n_dofs,P2UMF,'o-k');
loglog(n_dofs,P2CG_1em16,'d-k');
hold on;
loglog(n_dofs,P2CG_1em10,'^-k');   %,n_dofs,P2CG_1em4,'d--k'
% hold on;
% loglog(n_dofs,P2CG_1em4,'s-k');

xlabel('Number of DoFs');
ylabel('Absolute error');
xlim([1e0 1e+8]);
ylim([1e-16 1e+0]);

lgd = legend('Schur, 10^{-16}','Schur, 10^{-10}');  %'Eq. (9), UMFPACK',
lgd.FontSize = 12;
lgd.Location = 'northeast';
set(gca,'FontSize',18);

ax=gca;
ax.XTick = [1e0 1e2 1e4 1e6 1e8 1e10];
ax.YTick = [1e-20 1e-16 1e-12 1e-8 1e-4 1e0];

hold on;
% UMF_bound=5E-18*n_dofs.^1.0;
% loglog(n_dofs, UMF_bound,'--k');
% text(0.5*n_dofs(8),0.8*UMF_bound(6),'5\times10^{-18}\timesN^{1.0}','HorizontalAlignment','left','Color','k','FontSize',14);

CG_bound=2E-16*n_dofs.^1.0;
loglog(n_dofs, CG_bound,'--k');
text(0.13*n_dofs(1),1*CG_bound(7),'2\times10^{-16}\timesN^{1.0}','HorizontalAlignment','left','Color','k','FontSize',14);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'Pois_MM_UMFvsCG_grad','-dpdf','-r0');