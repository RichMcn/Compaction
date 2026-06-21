%PlotFig3bandc plots figures b and c. b shows the rescaled solution from a,
%showing collapse of the curves when Pcal is scaled by Scal. c shows
%collapse of solution from a when plotted with a solution using the same
%parameter except using Lcal/2 as the driving pressure. Also shows how
%changing Tcal has minimal effect on flow rate solutions

%Adds path for the solvers
addpath('Solvers')
z=linspace(0,1,1000);

%Parameters for baseline solution
Pcal = logspace(-3,5,1000);
Lcal =0; 
Scal = logspace(3,6,4);
Tcal=10^(-4);
Phi=0.05;

%Solves for baseline flow rates
[~,FlowRates] = VelLoopF(Pcal,Lcal,Scal,Tcal,Phi);

%Solves for gravity-driven solution
Pcal = 0;
Lcal = logspace(-3,5,1000);
[~,FlowRates2] = VelLoopF(Pcal,Lcal,Scal,Tcal,Phi);

%Solves for solution using a smaller value of Tcal
Pcal = logspace(-3,5,1000);
Lcal =0; 
Tcal=10^(-7);
[~,FlowRates3] = VelLoopF(Pcal,Lcal,Scal,Tcal,Phi);


% Define colors manually as an Nx3 RGB matrix (N = number of lines)
colors = [
    0    0.4470    0.7410;  % blue
    0.8500    0.3250    0.0980; % red-orange
    0.9290    0.6940    0.1250; % yellow
    0.4940    0.1840    0.5560; % purple
    0.4660    0.6740    0.1880; % green
];


subplot(2,2,[1,2])

%Plots scaled solution 3b
for i=1:4
    loglog(Pcal/Scal(i), -FlowRates(i,:), 'Color', colors(i,:), 'LineWidth', 1.5)
    hold on
end

%Formats plot
h = xlabel('$\mathcal{P}/\mathcal{S}$','interpreter','latex');
set(h, 'Units', 'normalized');              % Use normalized units for positioning
set(h, 'Position', [0.5, 1.15, 0])
text(1.5*10^(-8),2*10^(-4),"(b)", 'fontsize', 16) 
ylabel('Q_{z0}')
xlim([10^(-7),2*10^(-1)])
ylim([0.3*10^(-4),2*10^(-4)])
set(gca, 'fontsize', 16)
ax = gca;
ax.TickLength = [0.025 0.025];
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

%Plots comparison figure 3c
subplot(2,2,[3,4])

%Insible plots for legend
loglog([nan,nan],[nan,nan],'k', 'LineWidth', 1.5)
hold on
loglog([nan,nan],[nan,nan],'k--', 'LineWidth', 1.5)
loglog([nan,nan],[nan,nan],'k:', 'LineWidth', 1.5)


figure(1)

%Plots comparison curves
for i=1:length(Scal)

    loglog(Pcal, -FlowRates(i,:), 'Color', colors(i,:), 'LineWidth', 1)

    loglog(Lcal/2, -FlowRates2(i,:),'--', 'Color', colors(i,:), 'LineWidth', 2)

    loglog(Pcal, -FlowRates3(i,:),':', 'Color', colors(i,:), 'LineWidth', 1.5)

end

%Formats figure
box off
text(1.8*10^(-1),0.2*10^(-4),"(c)", 'fontsize', 16)
h = xlabel('$\mathcal{W}$','interpreter','latex');
set(h, 'Units', 'normalized');              % Use normalized units for positioning
set(h, 'Position', [0.5, 0.000000, 0])
ylabel('Q_{z0}')
ylim([0.2*10^(-4),0.2*10^(-3)])
xlim([10^(-1),10^5])
set(gca, 'fontsize', 16)
legendEntries = [{'$\mathcal{W} = \mathcal{P}$, $\mathcal{L}=0$, $\mathcal{T}=10^{-4}$'}, ...
                     {'$\mathcal{W} = \mathcal{L}/2$, $\mathcal{P}=0$, $\mathcal{T}=10^{-4}$'},...
                     {'$\mathcal{W} = \mathcal{P}$, $\mathcal{L}=0$, $\mathcal{T}=10^{-7}$'}];
hLeg = legend(legendEntries, 'Location', 'south');
set(hLeg, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 14);
legend boxoff
ax = gca;
ax.TickLength = [0.025 0.025];
hold off






