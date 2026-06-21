%PlotFig3a calculates and plots flow rate against pressure for Lcal=0 using
%VelLoopF. Also calculates radial profiles.

%Adds path for the solvers
addpath('Solvers')

%This figure is for Lcal=0
Pcal = logspace(-3,5,1000);
Lcal =0; 
Scal = logspace(3,6,4);
Tcal=10^(-4);
Phi=0.05;
z=linspace(0,1,1000);

%Solves for the flow rates
[~,FlowRates] = VelLoopF(Pcal,Lcal,Scal,Tcal,Phi);


% Define colors manually as an Nx3 RGB matrix (N = number of lines)
colors = [
    0    0.4470    0.7410;  % blue
    0.8500    0.3250    0.0980; % red-orange
    0.9290    0.6940    0.1250; % yellow
    0.4940    0.1840    0.5560; % purple
    0.4660    0.6740    0.1880; % green
];

%Calculates the radial profiles for inset
PSamples = [800,250,10];
for i=1:3
    [Chi,~,~] = MembraneProfApproxSol(PSamples(i),Lcal,Scal(2),Tcal,z);
        Fs(i,:)=Chi;
end
QSamples = interp1(Pcal,-FlowRates(2,:),PSamples);

%Plots main figure
figure(1)
loglog(Pcal, pi*Pcal*Phi^3/(1-Phi)^2,'k--','LineWidth',1.5)  % Darcy line: black dashed
hold on

%Hewitt plateau
loglog(Pcal, -ones(1,length(Pcal))*Phi*4*pi*(log(1-Phi)+Phi)/(9*(1-Phi)^2), 'Color', [0.8 0.8 0.8], 'LineStyle', '-.','LineWidth',2) % second line: grey dash-dot

%Invisible plot for the legend
loglog([10^(-9),10^(-8)],[100,100],'k:','LineWidth',1.5)

%Plots curves from model
for i=1:length(Scal)
    loglog(Pcal, -FlowRates(i,:), 'Color', colors(i,:), 'LineWidth', 1.5)
end

%Plots symbols for radial profiles
loglog(PSamples(1),QSamples(1),'ko', 'LineWidth', 1.5)
loglog(PSamples(2),QSamples(2),'k+', 'LineWidth', 1.5)
loglog(PSamples(3),QSamples(3),'kx', 'LineWidth', 1.5)

%Plots asymptotic solutions
for i=1:length(Scal)
    Qz0Asymp=nan(1,length(Pcal));
    for j=1:length(Pcal)
       [mxstr,Qz0Asymp(j)] = AsympFR(Phi,Phi./(1+9*Pcal(j)*(1-Phi)/4),Pcal(j),Scal(i), Tcal);
       if mxstr>0.2
           break;
       end
    end
   loglog(Pcal, Qz0Asymp,':', 'Color', colors(i,:), 'LineWidth', 1.5)
end
text(6*10^(-2),2*10^(-4),'$\pi \mathcal{P}\Phi^3/(1-\Phi)^2$', 'fontsize', 16,'interpreter','latex')
box off

%Formats figure
xlabel('$\mathcal{P}$','interpreter','latex')
ylabel('Q_{z0}')
ylim([7e-6, 2e-4])
yticks([10^(-5),10^(-4)])
xlim([1.7e-2, 1.1e5])
set(gca, 'fontsize', 16)
text(0.8*10^(-2),6*10^(-5),"(a)", 'fontsize', 16)
legendEntries = [{'Flow rate through stiff medium'}, ...
                     {'Plateau from Hewitt et al. (2016)'}, ...
                     {'Asymptotic solution'}];
DsStrings = arrayfun(@(D) sprintf('$\\mathcal{S} = 10^{%d}$', round(log10(D))), Scal, 'UniformOutput', false);
legendEntries = [legendEntries, DsStrings(1:end)];
hLeg = legend(legendEntries, 'Location', 'Best');
set(hLeg, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 16);
legend boxoff
ax = gca;
ax.TickLength = [0.025 0.025];
hold off

%Plots radial profiles in an inset figure
axes('Position',[0.685,0.32,0.2,0.2])
for i=1:3
    plot(Fs(i,:),z,'k','LineWidth',1)
    hold on
end
zpthalf = round(length(z)/2);
plot(Fs(1,zpthalf),z(zpthalf),'ko','LineWidth',1.5)
plot(Fs(2,zpthalf),z(zpthalf),'k+','LineWidth',1.5)
plot(Fs(3,zpthalf),z(zpthalf),'kx','LineWidth',1.5)
 xlim([-0.01,0.13])
 xlabel('$\chi(z)$','interpreter','latex')
ylabel('z')
set(gca, 'fontsize', 14)
ax = gca;
ax.TickLength = [0.05 0.05];
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';