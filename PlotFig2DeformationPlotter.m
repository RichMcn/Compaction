% PlotFig2DeformationPlotter solves for the deformation, porosity, and Eulerian
% Darcy velocity and outputs a Fig2 style plot for the choice of parameters
% entered in the first block

%Enter the model parameters to plot the deformation figure for that
%combination
Phi=0.05;
Pcal=1;
Lcal=1;
Scal=10;
Tcal=10^(-4);

%Plotting parameters
zpts=500;
epsilon=0.2;
rpts=1001;
stripes=10;

%Adds path for the solvers
addpath('Solvers')
lightGrey = [0.7, 0.7, 0.7];  % RGB for light grey

%Call solution
[maxstrain,Q,F,dFdz,phi0,dphi0dz,z] = NewtSolve_ApproxChi(Pcal,Lcal,Scal,Tcal,Phi,zpts);

%Initialise other parameters
dz = z(2)-z(1);
Dz0 = zeros(1,zpts);
r= linspace(0,1,rpts);
r_def=[];

%Deformed radial coordinates
for i=1:zpts
    r_def = [r_def;r*(1+F(i))];
end

%Displacements
for i=2:zpts
  Dz0(i) = trapz((phi0(1:i)'-Phi)./(1-phi0(1:i)')-2*F(1:i))*dz;
end

%Deformed vertical coordinates
z_def = z+Dz0';
z_defgrid=[];

%Vertical component of strain
eps_zzvec = (phi0'-Phi)./(1-phi0')-2*F;
eps_zz=[];
for i=1:rpts
    eps_zz = [eps_zz,eps_zzvec];
    z_defgrid = [z_defgrid,z_def];
end

%Eulerian Darcy velocity
[R,Z] =meshgrid(r,z);
ufield_r = Q * epsilon*R.*dFdz.*(1-phi0')./(1-Phi);
ufield_z = (Q * (Z./Z).*(1+(phi0'-Phi)/(1-Phi) -2*F)).*(1-phi0')./(1-Phi);
u_zEul=Q./((1+F).^2);


figure(1)
t = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact'); % Tight layout

% Large plot across multiple tiles
nexttile([2 3])
hold on
plot(epsilon*(1+F),z_def,'k');
plot(-epsilon*(1+F),z_def,'k');
plot([-epsilon,epsilon],[z_def(end),z_def(end)],'k')
xlim([-(2*epsilon)+0.05,2*epsilon+0.05;])
ylim([0,1])
set(gca, 'fontsize', 16)
xlabel('\epsilon r')
ylabel('z')
cb = colorbar('location','east');

title(cb, '$\varepsilon_{zz}$','interpreter','latex','fontsize',16);

contourf(epsilon*r_def,z_defgrid,eps_zz,25)
contourf(-epsilon*r_def,z_defgrid,eps_zz,25)

h = streamslice(epsilon*[-flip(r_def(:,2:end),2),r_def],[flip(z_defgrid(:,2:end),2),z_defgrid],[-flip(ufield_r(:,2:end),2),ufield_r],[flip(ufield_z(:,2:end),2),ufield_z],0.25);

set(h, 'Color', lightGrey, 'LineWidth', 1.5);

M=Scal;
delta = Tcal^(1/4);
alpha = (M+2)/M;
Omega = (alpha/4)^(1/4)/delta;
[phi0asymp] = Asympphi(Pcal,M,Phi,phi0(1),z,Omega);

% Smaller plot on the side
nexttile(4, [2 1])
plot(phi0,z_def,'LineWidth',1.5)
hold on
if Lcal==0
    plot(phi0asymp,z_def,'k:','LineWidth',1.5)
end
plot([Phi,Phi],[0,z_def(end)],'k--')
plot([phi0(1),phi0(1)],[0,z_def(end)],'k--')
yticks([])
ylim([0,1])
if Pcal>10
    xlim([-0.01,max(phi0)+0.01])
else
    xlim([0,max(phi0)+0.01])
end
set(gca, 'fontsize', 16)


xlabel('$\phi_0$','Interpreter','latex')
text(Phi-0.005,z_def(end)+0.05,'$\Phi$','Interpreter','latex','fontsize',16)
text(phi0(1)-0.005,z_def(end)+0.05,'$\varphi$','Interpreter','latex','fontsize',16)


% Smaller plot on the side
nexttile(5, [2 1])
plot(-u_zEul,z_def,'LineWidth',1.5)
yticks([])
ylim([0,1])
set(gca, 'fontsize', 16)
xlabel('$|u_z|$','interpreter','latex')
ax = gca;
ax.XAxisLocation = 'top';

% Ensure the x-axis label stays at the bottom
ax.XLabel.Position(2) = ax.YAxis.Limits(1) - 0.05 * diff(ax.YAxis.Limits);
ax.XLabel.VerticalAlignment = 'top';

