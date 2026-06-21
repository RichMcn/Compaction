function [Pcal,FlowRates] = VelLoopF(Pcal,Lcal,Scal,Tcal,Phi)
% VelLoopF
%
% This function calculates and plots the flow rates for a range of pressures
% and stiffness ratios
% The model stops calculation when the maximum strain exceeds 0.2
%
% The flow rates are computed using SolveUntilTol, and the results are plotted
% as log-log curves of pressure vs flow rate 


% Define pressure sweep 
    ptsMaterial = length(Scal);
     ptsPress = max(length(Pcal),length(Lcal));
if (isscalar(Pcal))&&(isscalar(Lcal))
    ptsPress = 300;
    Pcal = logspace(-3,1,ptsPress);
end

% Preallocate array to store computed flow rates
FlowRates = nan(ptsMaterial, ptsPress);  % Rows: D values, Cols: P values

% Loop over all combinations
for i = 1:ptsMaterial
    for j = 1:ptsPress

        % Solve the model for current parameter combination. Tolerance is
        % 1e-4
        if length(Pcal)>length(Lcal)
                 [maxstrain, Q, ~, ~, ~, ~, ~] = SolveUntilTol(Pcal(j), Lcal , Scal(i), Tcal, Phi, 1e-4, 1);
   
        else
                 [maxstrain, Q, ~, ~, ~, ~, ~] = SolveUntilTol(Pcal, Lcal(j), Scal(i), Tcal, Phi, 1e-4, 1);
           
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

end