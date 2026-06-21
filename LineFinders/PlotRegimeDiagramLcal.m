% Plots the regime curves for Lcal vs Scal

blue        = [0, 114, 178]   / 255;  % #0072B2
vermillion  = [213, 94, 0]    / 255;  % #D55E00
yellow      = [240, 228, 66]  / 255;  % #F0E442
green       = [0, 158, 115]   / 255;  % #009E73
lightGrey = [0.7, 0.7, 0.7];  % RGB for light grey

figure(1)
loglog(Scal(Scal>max(CompStrainCurve)),DarcyCurve(Scal>max(CompStrainCurve)),'Color', blue, 'LineWidth', 2)
hold on
loglog(Scal,StrCurves(5,:),'k')
loglog(Scal,StrCurves(1:4,:),'k:')
loglog(CompStrainCurve,Lcal,'k')
loglog(Scal(Scal>max(CompStrainCurve)),PlateauBottoms(1,(Scal>max(CompStrainCurve))),'Color', vermillion, 'LineWidth', 2)
loglog(PlateauTops(1,:),Lcal,'Color', vermillion, 'LineWidth', 2)
loglog(Scal(Scal>10^2),PlateauBottoms(2:3,Scal>10^2),'k', 'LineWidth', 1)
loglog(PlateauTops(2:3,Lcal>1),Lcal(Lcal>1),'k', 'LineWidth', 1)
loglog(Scal,compactioncurve,'k-.','LineWidth', 1.5)
loglog([10^(3),10^(4)],[10^(-2),10^(-2)],'Color', green, 'LineWidth', 2)
xlim([Scal(1),Scal(end)])
ylim([Lcal(1),Lcal(end)])
yticks([10^(-4),10^(-1),10^(2),10^5])
xticks([10^0,10^3,10^6])
xlabel("$\mathcal{S}$",'Interpreter','latex')
ylabel("$\mathcal{L}$",'Interpreter','latex')
set(gca,'fontsize',24)
set(gca,'TickDir','out');

colors = [
    0    0.4470    0.7410;  % blue
    0.8500    0.3250    0.0980; % red-orange
    0.9290    0.6940    0.1250; % yellow
    0.4940    0.1840    0.5560; % purple
    0.4660    0.6740    0.1880; % green
    213/255, 94/255, 0   
];



