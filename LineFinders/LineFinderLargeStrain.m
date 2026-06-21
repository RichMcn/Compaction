function [StressCurve] = LineFinderLargeStrain(Scal, Tcal, Pcal, Lcal, Phi, strainlimit)
% LineFinderLargeStrainNewt
%
% Identifies the boundary in parameter space where the maximum strain
% component |\varepsilon_{zz}| exceeds the specified strainlimit input.
%
% The method sweeps through stiffness ratios and a
% flow-driving parameter (Pcal or Lcal), checking at each point if the maximum
% strain exceeds the threshold.
%
% Inputs:
%   Scal          - Vector of stiffness ratios
%   Tcal          - Membrane bending parameter
%   Pcal          - Vector or scalar of pressure values
%   Lcal          - Vector or scalar of gravity parameter values
%   Phi         - Initial porosity
%   strainlimit - Threshold for |\varepsilon_{zz}|
%
% Outputs:
%   StressCurve  - Critical pressure or gravity parameter at which strain
%                    first exceeds 'strainlimit' for each material parameter

addpath('../Solvers')

% Lengths of input arrays
l_Ds = length(Scal);
l_Ps = length(Pcal);
l_Ls = length(Lcal);

% Sanity check: only one of each mechanical and flow-driving parameter can be swept
if  ((l_Ps > 1) && (l_Ls > 1))
    error('Only one of Ps or Ls, should be a vector.');
end

    T = Tcal;

% Set sweep direction: either Ps or Ls
if l_Ps > 1
    pressures = Pcal;
    L = Lcal;
else
    pressures = Lcal;
    P = Pcal;
end

% Initialize output index array
Pindices = nan(1, length(Scal));
jstart = 1;


% Sweep over mechanical parameter
for i = 1:length(Scal)

    % Stop if we're at the top of pressure range
    if jstart > length(pressures)-3
        break;
    end

    % Ensure index stays within bounds
    if jstart < 1
        jstart = 1;
    end

    % Sweep over flow parameter (pressure or gravity)
    for j = jstart:length(pressures)
     
            D = Scal(i);

        % Update pressure or gravity
        if l_Ps > 1
            P = Pcal(j);
        else
            L = Lcal(j);
        end

        % Call solver to compute max strain
        [maxstrain, ~, ~, ~, ~, ~, ~] = SolveUntilTol(P, L, D, T, Phi, 1e-4, 1);

        if isnan(maxstrain)
            break
        end


        % Check if strain exceeds threshold
        if maxstrain > strainlimit
            if j == 1
                break;  % Already exceeded at the first point
            end
            Pindices(i) = j - 1;
            jstart = j - 2;  % Speed up search by skipping ahead
            break;
        end

        % Optional progress display
        [i, j]

    end
end

% Build final stress curve from recorded indices
StressCurve = nan(size(Pindices));
valid = ~isnan(Pindices);
StressCurve(valid) = pressures(Pindices(valid));

end
