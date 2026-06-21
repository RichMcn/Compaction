% ConvergenceTester.m
% This script tests the convergence of a Newton solver for porosity and deformation
% equations by refining the spatial mesh and monitoring changes in the solution and flow-rate.

%Adds path for the solvers
addpath('Solvers')

%Parameter set under which convergence is tested. Convergence will be
%faster for Pcal+Lcal<1 as the boundary layer does not form
Pcal=10^(2);
Scal=10^4;
Tcal=10^(-4);
Phi=0.05;
Lcal=0;

%1=test approx solver, 0=test full solver
App=1;

%Initialised set of discretisation values for the convergence test
numDataPts=100;
minMesh=50;
maxMesh=100000;
zpts = round(logspace(log10(minMesh),log10(maxMesh),numDataPts));

%Initalises solution variables
Qvals = nan(numDataPts,1);
phisCompare=nan(numDataPts,maxMesh);

%Loops over the disretisation values calculating the solution each time
for i=1:length(zpts)

    if App==1
      [maxstrain,Q,F,dFdz,phi0,dphi0dz,z] = NewtSolve_ApproxChi(Pcal,Lcal,Scal,Tcal,Phi,zpts(i));
    else
          [maxstrain,Q,F,dFdz,phi0,dphi0dz,z] = NewtSolve(Pcal,Lcal,Scal,Tcal,Phi,zpts(i));
    end
     Qvals(i) = Q;
     phisCompare(i,:) = interp1(z,phi0,linspace(0,1,zpts(end)));

     %Prints to the screen current place in the loop
     PercentDone = i/length(zpts)
end

%Prints figure showing convergence of the infinity norm and 2-norm of the
%scheme for the solution phi_0 with the convergence of the solution for the
%flow rate Q plotted as an inset
figure(1)
loglog(zpts, vecnorm((phisCompare-phisCompare(end,:))',inf)/norm(phisCompare(end,:),inf )  )
hold on
loglog(zpts, vecnorm((phisCompare-phisCompare(end,:))',2)/norm(phisCompare(end,:),2 )  )
xlabel("M")
ylabel('||𝛷_M - 𝛷_{end}||_2 / ||𝛷_M||_2', 'Interpreter', 'tex')
set(gca,'fontsize', 14)
axes('Position',[0.5,0.5,0.3,0.3])
loglog(zpts,abs((Qvals-Qvals(end))/Qvals(end)))
xlabel("M")
ylabel("$(Q-Q_{end})/Q_{end}$",'interpreter', 'latex')
set(gca,'fontsize', 16)