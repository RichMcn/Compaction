function [pressurecurve] = LineFinderPlateauNewtBottom(Scal, Tcal, Pcal, Lcal, Phi, tol)
% LineFinderPlateauNewtBottom
% Identifies the transition line in parameter space where the flow rate approaches
% the compacted/plateau regime, i.e., becomes within 'tol' of the fully compacted limit.
%
% This occurs when stiffness ratio and flow drivers (P or L) are large.
%
% Inputs:
%   Scal       - Stiffness ratio
%   Tcal       - Membrane bending parameter
%   Pcal, Lcal     - Vector or scalar of pressure / gravity parameter
%   Phi        - Initial porosity
%   tol        - Tolerance for plateau approach (relative error)
%
% Outputs:
%   pressurecurve  - Value of Pcal or Lcal at which the flow enters the plateau regime

addpath('../Solvers')

% Parameter lengths
l_Ps = length(Pcal);
l_Ls = length(Lcal);

% Safety check: only one mechanical and one flow parameter should be swept
if (l_Ps > 1 && l_Ls > 1)
    error('Only one of Ds or Ts, and only one of Ps or Ls, should be a vector.');
end

% Compute target compacted flow rate (large parameters)
[~, targetQ, ~, ~, ~, ~, ~] = SolveUntilTol(1e3, 1e3, 1e12, 1e12, Phi, 1e-4, 1);



% Identify pressure/gravity sweep
if l_Ps > 1
    pressures = Pcal;
    L = Lcal;
else
    pressures = Lcal;
    P = Pcal;
end

% Start search above P > 1 (below which compaction doesn't occur)
jstart = find(pressures > 1, 1);
if isempty(jstart)
    pressurecurve = nan(1, length(Scal));
    return
end

Pindices = nan(1, length(Scal));
search_done = false;

% Sweep through material parameters
for i = length(Scal):-1:1
    if jstart < 1 || jstart >= length(pressures)
        break
    end

    j = jstart;

    while j < length(pressures)


            D = Scal(i);

        % Assign P or L
        if l_Ps > 1
            P = Pcal(j);
        else
            L = Lcal(j);
        end

        % Compute flow rate
        [maxstrain, Q, ~, ~, ~, ~, ~] = SolveUntilTol(P, L, D, Tcal, Phi, 1e-4, 1);

        if isnan(Q)
            search_done=1;
            break
        end

        rel_err = abs((targetQ - Q) / Q);

        % If we've already reached the plateau at the start, go further down
        if j == jstart && rel_err < tol
            j = j - 2;
            jstart = max(j - 2, 1);
            continue
        end

        % Record the last value before hitting plateau
        if rel_err < tol
            Pindices(i) = j - 1;
            jstart = max(j - 2, 1);
            break
        end

        % Exit if strain grows too large (invalid model regime)
        if maxstrain > 0.2
            search_done = true;
            break
        end

        disp([i, j])
        j = j + 1;
    end

    if search_done
        break
    end
end

% Build pressure curve
pressurecurve = nan(size(Pindices));
valid = ~isnan(Pindices);
pressurecurve(valid) = pressures(Pindices(valid));

end
