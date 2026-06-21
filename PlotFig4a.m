%PlotFig4a plots flow rate against Pcal figures and radial profiles for smaller values
%of Scal for figure 4a with Lcal=0

% Define colors manually as an Nx3 RGB matrix (N = number of lines)
colors = [
    0    0.4470    0.7410;  % blue
    0.8500    0.3250    0.0980; % red-orange
    0.9290    0.6940    0.1250; % yellow
    0.4940    0.1840    0.5560; % purple
    0.4660    0.6740    0.1880; % green
];

%Adds path for the solvers
addpath('Solvers')

%This figure is for Lcal=0
Pcal = logspace(-3,2,1000);
Lcal =0; 
Scal = [0.2,1,10,100];
Tcal=10^(-4);
Phi=0.05;

%Calculates flow rates
[~,FlowRates] = VelLoopF(Pcal,Lcal,Scal,Tcal,Phi);

%Calculates radial profiles
z = linspace(0,1,1000);
PSamples = [ 0.2,0.1,0.05;1,0.2,0.1 ];
for i=1:3
    [Chi,~,~] = MembraneProfApproxSol(PSamples(1,i),Lcal,Scal(1),Tcal,z);
        Fs1(:,i)=Chi;
end
for i=1:3
    [Chi,~,~] = MembraneProfApproxSol(PSamples(2,i),Lcal,Scal(3),Tcal,z);
        Fs2(:,i)=Chi;
end
QSamples = [interp1(Pcal,-FlowRates(1,:),PSamples(1,:));interp1(Pcal,-FlowRates(3,:),PSamples(2,:)) ];


%Plots main flow rate figure
figure(1)
loglog(Pcal,-FlowRates,'LineWidth',1.5)

%Formates figure
 ylim([0.3*10^(-5),2*10^(-4)])
 xlim([0.6*10^(-2),0.3*10^2])
hold on
box off
loglog(Pcal, pi*Pcal*Phi^3/(1-Phi)^2,'k--','LineWidth',1.5)  % first line: black dashed
text(2.5*10^(-1),2*10^(-4),'\pi P\Phi^3/(1-\Phi)^2', 'fontsize', 16)
xlabel('$\mathcal{P}$','Interpreter','latex')
ylabel('Q_{z0}')
set(gca, 'fontsize', 16)
ScalStrings = arrayfun(@(D) sprintf('$\\mathcal{S} = 10^{%d}$',...
    round(log10(D))), Scal, 'UniformOutput', false);
legendEntries = [ScalStrings(1:end)];
hLeg = legend(legendEntries, 'Location', 'Northwest');
set(hLeg, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 16);
xticks([10^(-2),1])
legend boxoff
ax = gca;
ax.TickLength = [0.025 0.025];

%Plots symbols to link to inset
loglog(PSamples(1,1),QSamples(1,1),'ko', 'LineWidth', 1.5,'HandleVisibility','off')
loglog(PSamples(1,2),QSamples(1,2),'k+', 'LineWidth', 1.5,'HandleVisibility','off')
loglog(PSamples(1,3),QSamples(1,3),'kx', 'LineWidth', 1.5,'HandleVisibility','off')

loglog(PSamples(2,1),QSamples(2,1),'k*', 'LineWidth', 1.5,'HandleVisibility','off')
loglog(PSamples(2,2),QSamples(2,2),'kdiamond', 'LineWidth', 1.5,'HandleVisibility','off')
loglog(PSamples(2,3),QSamples(2,3),'kpentagram', 'LineWidth', 1.5,'HandleVisibility','off')

hold off

%Plots first inset
axes('Position',[0.52,0.32,0.2,0.2])
for i=1:3
    plot(Fs1(:,i),z,'Color', colors(1,:),'LineWidth',1.5)
    hold on
end
zpthalf = round(length(z)/2);
plot(Fs1(zpthalf,1),z(zpthalf),'ko','LineWidth',1.5)
plot(Fs1(zpthalf,2),z(zpthalf),'k+','LineWidth',1.5)
plot(Fs1(zpthalf,3),z(zpthalf),'kx','LineWidth',1.5)

%Formats
 xlim([-0.01,0.12])
 xlabel('\chi(z)')
ylabel('z')
set(gca, 'fontsize', 14)
ax = gca;
ax.TickLength = [0.05 0.05];
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

%Plots second inset
axes('Position',[0.785,0.32,0.2,0.2])
for i=1:3
    plot(Fs2(:,i),z,'Color', colors(3,:),'LineWidth',1.5)
    hold on
end
zpthalf = round(length(z)/2);
plot(Fs2(zpthalf,1),z(zpthalf),'k*','LineWidth',1.5)
plot(Fs2(zpthalf,2),z(zpthalf),'kdiamond','LineWidth',1.5)
plot(Fs2(zpthalf,3),z(zpthalf),'kpentagram','LineWidth',1.5)

%Formats
 xlim([-0.005,0.13])
 xlabel('\chi(z)')
set(gca, 'fontsize', 14)
ax = gca;
ax.TickLength = [0.05 0.05];
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';