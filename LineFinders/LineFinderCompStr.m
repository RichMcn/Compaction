function [ScalCurve] = LineFinderCompStr(Scal, Tcal, Pcal, Lcal, Phi)
% LineFinderCompStr
% ----------------------------------------------
% This function identifies a curve in parameter space for the regime map
% which shows where the contact condition fails. The line is close to
% vertical showing that this condition is mainly a function of Scal
%
% Inputs:
%   Scal    - Vector of stiffness ratios
%   Tcal    - Membrane bending parameter
%   Pcal    - Vector of pressure values
%   Lcal    - Vector of L values (one of Pcal or Lcal should be a single
%   value
%   Phi   - Initial porosity
%   tol   - Tolerance for the relative error between 
%
% Output:
%   ScalCurve - Curve of Scal values (to be plotted against Pcal or Lcal)
%

addpath('../Solvers')

% --- Determine input dimensions
l_Ds = length(Scal);
l_Ps = length(Pcal);
l_Ls = length(Lcal);

% --- Check input sweep compatibility
if ((l_Ps > 1) && (l_Ls > 1))
    return
end

% --- Set up material parameters and geometry
    materialparams = Scal;
    T = Tcal;

if l_Ps > 1
    pressures = Pcal;
    L = Lcal;
else
    pressures = Lcal;
    P = Pcal;
end

% --- Initialise index array
Scalindices = nan(1, length(pressures));



% --- Flags
finish = false;

jstart=1;

% --- Main loop over pressure values
for i = 1:length(pressures)


    if jstart >= (length(materialparams) - 3)
        break;
    end

    j = max(jstart, 1);
    dropfurther = false;

    while j < length(materialparams)

              [i,j]

        % Assign current D and T
        if l_Ds > 1
            D = Scal(j);
        else
            T = Tcal(j);
            D = 1e-2 * Tcal(j); % Empirical scaling
        end

        % Assign current P and L
        if l_Ps > 1
            P = Pcal(i);
        else
            L = Lcal(i);
        end

        % Solve for Q and strain

        [maxstrain,~,F,~,phi0,~,z] = SolveUntilTol(P,L,D,T,Phi,1e-4,1);

        % Stop if max strain exceeds limit
        if maxstrain > 0.2
            finish = true;
        
            break;
        end

        MinCompStr = min(P+L*(1-z)+(2/3)*(phi0'-Phi)/(1-Phi) - 2*F);

 

        if (MinCompStr > 0) && (j == jstart)
            j = max(j - 5, 1);
            jstart = max(j - 5, 1);
            continue;
        end

        if MinCompStr > 0
            if j < 2
                break;
            end

            if dropfurther
                j = max(j - 5, 1);
                jstart = max(j - 5, 1);
                continue;
            end

            Scalindices(i) = j - 1;
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
valid = ~isnan(Scalindices);
ScalCurve = nan(1, length(Scalindices));
ScalCurve(valid) = materialparams(Scalindices(valid));


end
