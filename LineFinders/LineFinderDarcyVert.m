function [DarcyCurve] = LineFinderDarcyVert(Scal, Tcal, Pcal, Lcal, Phi,tol)
% LineFinderDarcyVertNewt
% ----------------------------------------------
% This function identifies a curve in parameter space
% where the flow rate Q falls below the analytical Darcy velocity. This
% finds a curve to be plotted versus the pressure parameter
%
% Inputs:
%   Scal    - Vector of stiffness ratios
%   Tcal    - Membrane bending parameter
%   Pcal    - Vector of pressure values
%   Lcal    - Vector of L values
%   Phi   - Porosity
%   tol   - Tolerance for the relative error between calculated flow rate
%   and Darcy flow rate
%
% Output:
%   DarcyCurve - Curve of Scal values to be plotted vs pressure
%

addpath('../Solvers')

% --- Determine input dimensions
l_Ps = length(Pcal);
l_Ls = length(Lcal);

% --- Check input sweep compatibility
if  ((l_Ps > 1) && (l_Ls > 1))
    return
end


if l_Ps > 1
    pressures = Pcal;
    L = Lcal;
else
    pressures = Lcal;
    P = Pcal;
end

% --- Initialise index array
Dindices = nan(1, length(pressures));

% --- Initial guess for jstart: first value above D ~ 1e-2
jstart = find(Scal > 1e-1, 1);

% --- Initial guess for istart: first pressure above ~1e-2
istart = find(pressures > 1e-2, 1);

% --- Compute the analytical Darcy flow rate
Qdarcy = Phi^3 * pi * pressures / (1 - Phi)^2;

% --- Flags
finish = false;

jstart=1;
j=1;

% --- Main loop over pressure values
for i = istart:length(pressures)


    if jstart >= (length(Scal) - 3)
        break;
    end

    j = max(jstart, 1);
    dropfurther = false;

    while j < length(Scal)

              [i,j]

      
            D = Scal(j);

        % Assign current P and L
        if l_Ps > 1
            P = Pcal(i);
        else
            L = Lcal(i);
        end

        % Solve for Q and strain
        [maxstrain, Q, ~, ~, ~, ~, ~] = SolveUntilTol(P, L, D, Tcal, Phi, 1e-4, 1);

        % Stop if max strain exceeds limit
        if maxstrain > 0.2
            finish = true;
            break;
        end

        % Check whether Darcy curve has been crossed
        rel_diff = ( 1 - (Qdarcy(i)/(-Q)) );

        if (rel_diff < tol) && (j == jstart)
            j = max(j - 5, 1);
            jstart = max(j - 5, 1);
            continue;
        end

        if rel_diff < tol
            if j < 2
                break;
            end

            if dropfurther
                j = max(j - 5, 1);
                jstart = max(j - 5, 1);
                continue;
            end

            Dindices(i) = j - 1;
            j=j-5;
            jstart = j - 5;
            dropfurther = true;
            break;
        end

      
        dropfurther = false;
        j = j + 1;
    end

    if finish
        break;
    end
end

% --- Convert valid indices to corresponding material parameter values
valid = ~isnan(Dindices);
DarcyCurve = nan(1, length(Dindices));
DarcyCurve(valid) = Scal(Dindices(valid));


end
