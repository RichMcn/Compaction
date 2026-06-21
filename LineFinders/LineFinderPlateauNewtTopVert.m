function [MaterialParamCurve] = LineFinderPlateauNewtTopVert(Scal, Tcal, Pcal, Lcal, Phi, tol)
% LineFinderPlateauNewtTopVert
% Finds the line in parameter space where the flow rate approaches the
% compacted plateau value from the top side (descending over pressures),
% for varying stiffness ratio.
%
% Inputs:
%   Scal       - Stiffness ratio
%   Tcal       - Bending parameter
%   Pcal, Lcal - vectors or scalars for pressure and gravity parameters
%   Phi        - initial porosity
%   tol        - tolerance for closeness to plateau flow rate
%
% Outputs:
%   MaterialParamCurve  - vector of material parameter values at which flow reaches plateau


addpath('../Solvers')

l_Ds = length(Scal);
l_Ps = length(Pcal);
l_Ls = length(Lcal);

% Check input validity: only one of Ds or Ts should be vector,
% and only one of Ps or Ls should be vector
if ((l_Ps > 1) && (l_Ls > 1))
    return
end

if l_Ps > 1
    pressures = Pcal;
    L = Lcal;
else
    pressures = Lcal;
    P = Pcal;
end

jstart = 1;

% Compute the target compacted flow rate for large parameters
[~, target, ~, ~, ~, ~, ~] = SolveUntilTol(1e3, 1e3, 1e12, 1e12, Phi, 1e-4, 1);

finish = false;

% Find starting indices where to begin descending through pressures
Stoploop = false;
for i = length(pressures):-1:1
    if Stoploop
        break;
    end
    for j = 1:length(Scal) - 1

        
            if Scal(j) > 100*pressures(i)
                Stoploop = true;
                istart = i;
                jstart = j;
                break
            end

    end
end

Dindices1 = nan(1, length(pressures));

% Main loop: sweep backwards over pressures starting from istart
for i = istart:-1:1
    if jstart < 1 || finish
        break;
    end
    
    j = jstart;
    
    while j < length(Scal)

         [i, j]

        % Early cutoff condition:
        % If near bottom of materialparam range and pressure is low,
        % we are beyond plateau region - stop further searching.
        if (j == length(Scal) - 2) && (pressures(i) < 1)
            finish = true;
            break;
        end
        
        if j < 1
            break;
        end
        
            D = Scal(j);
        
        if l_Ps > 1
            P = Pcal(i);
        else
            L = Lcal(i);
        end
        
        % Solve for strain and flow at current parameters
        [maxstrain, Q, ~, ~, ~, ~, ~] = SolveUntilTol(P, L, D, Tcal, Phi, 1e-4, 1);

        if isnan(Q)
            finish=1;
            break
        end
        
        % If we're close to the target plateau flow rate at the starting
        % j index, drop back two steps and try again to find exact boundary
        if (abs((Q - target) / Q) < tol) && (j == jstart)
            j = j - 2;
            jstart = j - 2;
            continue;
        end
        
        % If close enough to plateau and not the first index, record position
        if (abs((Q - target) / Q) < tol) && (j ~= jstart)
            Dindices1(i) = j - 1;
            jstart = j - 2;
                if maxstrain>0.2
                    finish;
                end
            break;
        end
        
    %    Uncomment if you want to break on strain threshold:
        if maxstrain > 0.2
            finish=1;
            break;
        end

        if j==length(Scal)-1
            jstart=length(Scal)-3;
        end
        
        % Optional debug printout to monitor progress:
       
        
        j = j + 1;
    end
end

% Construct output curve vector from indices
MaterialParamCurve = nan(size(Dindices1));
validIdx = ~isnan(Dindices1);
MaterialParamCurve(validIdx) = Scal(Dindices1(validIdx));


end
