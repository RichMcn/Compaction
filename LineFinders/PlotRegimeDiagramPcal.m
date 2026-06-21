% Plots the regime curves for Pcal vs Scal


blue        = [0, 114, 178]   / 255;  % #0072B2
vermillion  = [213, 94, 0]    / 255;  % #D55E00
yellow      = [240, 228, 66]  / 255;  % #F0E442
green       = [0, 158, 115]   / 255;  % #009E73
lightGrey = [0.7, 0.7, 0.7];  % RGB for light grey

figure(1)
loglog(Scal,compactioncurve,'-.','Color',lightGrey,'LineWidth', 1.5)
hold on
loglog(Scal(Scal>max(CompStrainCurve)),StrCurves(1:4,Scal>max(CompStrainCurve)),'k:', 'LineWidth', 1.5)
loglog(Scal(Scal>max(CompStrainCurve)),DarcyCurve(Scal>max(CompStrainCurve)),'Color', blue, 'LineWidth', 2)
loglog(DarcyVertCurve,Pcal,'Color', blue, 'LineWidth', 2)
loglog(DarcyVertCurve2,Pcal,'Color', blue, 'LineWidth', 2)
loglog(Scal,StrCurves(5,:),'k')
loglog(CompStrainCurve,Pcal,'k')
loglog(Scal(Scal>min(CompStrainCurve)),PlateauBottoms(1,(Scal>min(CompStrainCurve))),'Color', vermillion, 'LineWidth', 2)
loglog(PlateauTops(1,:),Pcal,'Color', vermillion, 'LineWidth', 2)
loglog(Scal(Scal>10^2),PlateauBottoms(2:3,Scal>10^2),'k', 'LineWidth', 1)
loglog(PlateauTops(2:3,Pcal>1),Pcal(Pcal>1),'k', 'LineWidth', 1)
loglog([2*10^(-1),5*10^(-1)],[10^(-1),10^(-1)],'Color', yellow, 'LineWidth', 2)
loglog([10^(3),10^(4)],[10^(-2),10^(-2)],'Color', green, 'LineWidth', 2)
xlim([5*10^(-2),Scal(end)])
ylim([Pcal(1),Pcal(end)])
xlabel("$\mathcal{S}$",'Interpreter','latex')
ylabel("$\mathcal{P}$",'Interpreter','latex')
set(gca,'fontsize',16)
set(gca,'TickDir','out');



colors = [
    0    0.4470    0.7410;  % blue
    0.8500    0.3250    0.0980; % red-orange
    0.9290    0.6940    0.1250; % yellow
    0.4940    0.1840    0.5560; % purple
    0.4660    0.6740    0.1880; % green
    213/255, 94/255, 0   
];




% The block below can be used to plot the inset figure and lay lines on the
% regime figure


% for i=1:5
% 
% [~,I]=max(-FlowRates(i,:));
% 
% loglog([ScalforCurves(i),ScalforCurves(i)],[PcalforCurves(1),PcalforCurves(I)] ,'Color', colors(i,:),'Linestyle','--' ,'LineWidth', 2 )
% 
% end
% 
% figure(2)
% 
% loglog(PcalforCurves,pi*PcalforCurves*Phi^3/(1-Phi)^2,'k--','LineWidth', 1.5 )
% hold on
% 
% for i=1:5
%     loglog(PcalforCurves,-FlowRates(i,:) ,'Color', colors(i,:) ,'LineWidth', 2 )
% end
% 
% xlim([0.3*10^(-1),0.3*10^6])
% ylim([0.3*10^(-4),2*10^(-4)])
% xlabel("$\mathcal{P}$",'Interpreter','latex')
% ylabel("Q_{z0}")
% set(gca,'fontsize',24)
% 
% text(6*10^(-2),6*10^(-4),'\pi P\Phi^3/(1-\Phi)^2', 'fontsize', 24)
% 
% ScalStrings = arrayfun(@(D) sprintf('$\\mathcal{S} = 10^{%d}$', round(log10(D))), ScalforCurves, 'UniformOutput', false);
% 
% legendEntries = ['Darcy',ScalStrings(1:end)];
% 
% hLeg = legend(legendEntries, 'Location', 'Best');
% set(hLeg, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 24);