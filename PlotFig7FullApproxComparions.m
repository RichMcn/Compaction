%PlotFig7FullApproxComparisons Plots a figure which outputs a plot showing
%the relative error in flow-rate between the solution of the full coupled
%system and the solution found using the approximation for Chi

%Adds path for the solvers
addpath('Solvers')

%Parameter conbinations over which we sweep
Pcal = logspace(-4,6,200);
ScalTcalPhi = [0.2,10^(-4),0.01;
                0.2,10^(-4),0.05;
                0.2,10^(-4),0.1;
                0.2,10^(-4),0.2;
                0.2,10^(-8),0.01;
                0.2,10^(-8),0.05;
                0.2,10^(-8),0.1;
                0.2,10^(-8),0.2;
                10,10^(-4),0.01;
                10,10^(-4),0.05;
                10,10^(-4),0.1;
                10,10^(-4),0.2;
                10^4,10^(-4),0.01;
                10^4,10^(-4),0.05;
                10^4,10^(-4),0.1;
                10^4,10^(-4),0.2];


[a,b] = size(ScalTcalPhi);
FRapprox = nan(a+4,length(Pcal));
FRfull = nan(a+4,length(Pcal));

%Calculates approximate solution
for i=1:a

    for j=1:length(Pcal)

        [1,i,j]

        [maxstrain,Q,F,dFdz,phi0,dphi0dz,z] = SolveUntilTol(Pcal(j),0,ScalTcalPhi(i,1),ScalTcalPhi(i,2),ScalTcalPhi(i,3),1e-3,1);
    
       if ScalTcalPhi(i,2)==10^(-4)
            if maxstrain>0.2
                break;
            end
       else
            if maxstrain>0.15
                break;
            end
       end

        FRapprox(i,j)=Q;

    end

end

%Also calculates some solutions using Lcal instead of Pcal
for i=1:4

        for j=1:length(Pcal)

            [1,i,j]
    
            [maxstrain,Q,F,dFdz,phi0,dphi0dz,z] = SolveUntilTol(0,Pcal(j),ScalTcalPhi(i,1),ScalTcalPhi(i,2),ScalTcalPhi(i,3),1e-3,1);
        
           if ScalTcalPhi(i,2)==10^(-4)
                if maxstrain>0.2
                    break;
                end
           else
                if maxstrain>0.15
                    break;
                end
           end
    
            FRapprox(a+i,j)=Q;

        end


end

%Calculates solutions from the full coupled system using same paramters as
%above
for i=1:a
  
    for j=1:length(Pcal)

           [2,i,j]

        [maxstrain,Q,F,dFdz,phi0,dphi0dz,z] = SolveUntilTol(Pcal(j),0,ScalTcalPhi(i,1),ScalTcalPhi(i,2),ScalTcalPhi(i,3),1e-3,0);

        if maxstrain>0.25
            break;
        end

        FRfull(i,j)=Q;

    end

end

%Lcal solutions
for i=1:4
  
    for j=1:length(Pcal)

           [2,i,j]

        [maxstrain,Q,F,dFdz,phi0,dphi0dz,z] = SolveUntilTol(0,Pcal(j),ScalTcalPhi(i,1),ScalTcalPhi(i,2),ScalTcalPhi(i,3),1e-3,0);

        if maxstrain>0.25
            break;
        end

        FRfull(a+i,j)=Q;

    end

end

%Plots comparison figure
Lstyles = ["-","--",":","-."];

for i=1:4
    if i==1
        loglog(Pcal,abs((FRapprox(i,:)-FRfull(i,:))./FRfull(i,:)),'b','linestyle',Lstyles(i),'LineWidth',1)
        hold on
    else
        loglog(Pcal,abs((FRapprox(i,:)-FRfull(i,:))./FRfull(i,:)),'b','linestyle',Lstyles(i),'LineWidth',1,'HandleVisibility','off')
    end
end

for i=5:8
    if i==5
        loglog(Pcal,abs((FRapprox(i,:)-FRfull(i,:))./FRfull(i,:)),'r','linestyle',Lstyles(i-4),'LineWidth',1)
    else
           loglog(Pcal,abs((FRapprox(i,:)-FRfull(i,:))./FRfull(i,:)),'r','linestyle',Lstyles(i-4),'LineWidth',1,'HandleVisibility','off')
    end
end

%Plot for legend
loglog([10^(-2),10^2],[100,100],'k-','LineWidth',1)

for i=13:16
    if i==13
        loglog(Pcal,abs((FRapprox(i,:)-FRfull(i,:))./FRfull(i,:)),'m','linestyle',Lstyles(i-12),'LineWidth',1)
    else
            loglog(Pcal,abs((FRapprox(i,:)-FRfull(i,:))./FRfull(i,:)),'m','linestyle',Lstyles(i-12),'LineWidth',1,'HandleVisibility','off')
    end
end

for i=17:20
    if i==17
        loglog(Pcal,abs((FRapprox(i,:)-FRfull(i,:))./FRfull(i,:)),'g','linestyle',Lstyles(i-16),'LineWidth',1)
    else
        loglog(Pcal,abs((FRapprox(i,:)-FRfull(i,:))./FRfull(i,:)),'g','linestyle',Lstyles(i-16),'LineWidth',1,'HandleVisibility','off')
    end
end

for i=9:12
loglog(Pcal,abs((FRapprox(i,:)-FRfull(i,:))./FRfull(i,:)),'k','linestyle',Lstyles(i-8),'LineWidth',1)
hold on
end

maxerr = max(abs((FRapprox-FRfull)./FRfull),[],'all');

loglog([10^(-4),10^4],[maxerr,maxerr],'k:')

%Formatting
 xlabel("$\mathcal{P}$",'Interpreter','latex')
 ylabel("Rel Err(Q,Q_{full})")
 ylim([10^(-10),1])
 xlim([10^(-4),10^5])
 xticks([10^(-3),1,10^3])
legend("$\mathcal{S}=0.2$, $\mathcal{T}=10^{-4}$","$\mathcal{S}=0.2$, $\mathcal{T}=10^{-8}$",...
    "$\mathcal{S}=10$, $\mathcal{T}=10^{-4}$","$\mathcal{S}=10^4$, $\mathcal{T}=10^{-4}$",...
    "$\mathcal{S}=0.2$, $\mathcal{T}=10^{-4}$ ($\mathcal{L}$)",...
    "$\Phi=0.01$","$\Phi=0.05$","$\Phi=0.1$","$\Phi=0.2$",'interpreter','latex','NumColumns',2)

 set(gca, 'fontsize', 16)
