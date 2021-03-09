clear all;
close all;
data = dlmread('data_error_grad.txt','',1,0);
n_dofs = data(:,1);
UMF = data(:,2);
M_1em16 = data(:,3);
M_1em10 = data(:,4);
M_UMF = data(:,5);

h = figure;
loglog(n_dofs,UMF,'o-k');
hold on;
loglog(n_dofs,M_1em16,'d-k');
loglog(n_dofs,M_1em10,'^-k');  % ,'Color',[0.4,0.4,0.4] 
loglog(n_dofs,M_UMF,'s-k');     %

xlabel('Number of DoFs');
ylabel('Absolute error');
xlim([1e0 1e+8]);
ylim([1e-16 1e+0]);

lgd = legend('Monolithic','{\ittol}_{\itprm}=10^{-16}','{\ittol}_{\itprm}=10^{-10}','UMFPACK');  
lgd.FontSize = 16;
lgd.Location = 'northeast';
set(gca,'FontSize',20);

ax=gca;
ax.XTick = [1e0 1e2 1e4 1e6 1e8 1e10];
ax.YTick = [1e-20 1e-16 1e-12 1e-8 1e-4 1e0];

hold on;
UMF_bound=5E-18*n_dofs.^1.0;
loglog(n_dofs, UMF_bound,'--k');
text(0.5*n_dofs(8),0.8*UMF_bound(6),'\alpha_R=5\times10^{-18}','HorizontalAlignment','left','Color','k','FontSize',18);

CG_bound=2E-16*n_dofs.^1.0;
loglog(n_dofs, CG_bound,'--k');
text(0.10*n_dofs(1),0.6*CG_bound(2),'\alpha_R=2\times10^{-16}','HorizontalAlignment','left','Color','k','FontSize',18,'Rotation',20);

tria_ref=n_dofs;
tria_bound=UMF_bound/100;

tria_start=16; tria_delta=3;
line([tria_ref(tria_start),tria_ref(tria_start+tria_delta)], [tria_bound(tria_start),tria_bound(tria_start)], 'Color', 'k');
line([tria_ref(tria_start+tria_delta),tria_ref(tria_start+tria_delta)], [tria_bound(tria_start),tria_bound(tria_start+tria_delta)], 'Color', 'k');
line([tria_ref(tria_start),tria_ref(tria_start+tria_delta)], [tria_bound(tria_start),tria_bound(tria_start+tria_delta)], 'Color', 'k');
text(1.2*tria_ref(tria_start+rem(tria_delta,2)),0.25*tria_bound(tria_start),'1','FontSize',18)
text(1.2*tria_ref(tria_start+tria_delta),1.2*tria_bound(tria_start+rem(tria_delta,2)),'1','FontSize',18)

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'Pois_UMFvsCG_MM_Schur_10tom16_other_variant_grad','-dpdf','-r0');