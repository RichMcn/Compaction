% VelLoop
%
% This script calculates and plots the flow rate for a range of pressures 
% and stiffness ratios
% The model uses a porosity value and stops calculation 
% when the maximum strain exceeds 0.2 (can be changed)
% The flow rates are computed using SolveUntilTol, and the results are plotted
% as log-log curves of pressure vs flow rate for different values of Scal

%Adds path for the solvers
addpath('Solvers')

% Define pressure sweep (logarithmic scale)
ptsPress = 1500;
ptsMaterial = 5;
Pcal = logspace(-4,5,ptsPress);

% Define membrane bending
Tcal =   1e-07;   

% Other model parameters
Scal = logspace(-2,2,ptsMaterial);   
Phi = 0.05;    % Porosity
Lcal = 0;     % Gravity parameter 

% Preallocate array to store computed flow rates
FlowRates = nan(ptsMaterial, ptsPress);  % Rows: D values, Cols: P values

% Loop over all combinations 
for i = 1:ptsMaterial
    for j = 1:ptsPress

        % Solve the model for current parameter combination
        if length(Pcal)>length(Lcal)
                 [maxstrain, Q, F, dFdz, phi0, dphi0dz, z] = SolveUntilTol(Pcal(j), Lcal , Scal(i), Tcal, Phi, 1e-4, 1);
          
        else
                 [maxstrain, Q, F, dFdz, phi0, dphi0dz, z] = SolveUntilTol(Pcal, Lcal(j), Scal(i), Tcal, Phi, 1e-4, 1);
        
        end

        % If the maximum strain exceeds threshold, stop evaluating this D
        if maxstrain > 0.2
            break;
        end

        % Store the flow rate for valid strain regime
        FlowRates(i, j) = Q;

        % Optional: display current indices (for debugging/tracking)
        [i, j]
    end
end

% Plot the results: log-log plot of pressure vs flow rate
figure(1)
loglog(Pcal, -FlowRates) 
hold on
loglog(Pcal,pi*Pcal*Phi^3/(1-Phi)^2,'k--')
xlabel('P')
ylabel('Q_{z0}')
ylim([10^(-6),10^(-3)])
xlim([10^(-3),10^5])
set(gca,'fontsize',16)
legend(arrayfun(@(D) sprintf('$\\mathcal{S} = %.0f$', D), Scal, 'UniformOutput', false), 'Location', 'Best','Interpreter','latex')