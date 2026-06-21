function [maxstrain,Q,Chi,dChidz,phi0,dphi0dz,z] = SolveUntilTol(Pcal,Lcal,Scal,Tcal,Phi,tol,Apr)
% SolveUntilTol solves the coupled membrane deformation and porosity equations 
% using a mesh refinement scheme until convergence is achieved within a specified tolerance.
%
% Inputs:
%   Pcal   - Applied pressure
%   Lcal   - Gravitational pressure parameter
%   Scal   - Ratio of stiffnesses
%   Tcal   - Membrane bending parameter
%   Phi - Initial porosity of the undeformed medium
%   tol - Relative convergence tolerance 
%   Apr - Flag (1 = use approximate solution for membrane profile, 0 = solve full coupled system)
%
% Outputs:
%   maxstrain  - Maximum strain in the porosity field
%   Q          - Volumetric flux
%   Chi          - Membrane deformation profile
%   dChidz       - Derivative of membrane deformation
%   phi0       - Porosity profile
%   dphi0dz    - Derivative of porosity
%   z          - Spatial mesh




%Gridpoints on which to test the mesh, doubling the size each time
zpts = [625,1250,2500,5000,10000,20000,40000];
phi0prev=Phi*ones(zpts(1),1);
Converged=0;

    % Loop over increasing mesh sizes to check for convergence
    for i=1:length(zpts)
        
        %Selects either the full system to solve or the approximate system
        if Apr
            [maxstrain,Q,Chi,dChidz,phi0,dphi0dz,z] = NewtSolve_ApproxChi(Pcal,Lcal,Scal,Tcal,Phi,zpts(i));
        else
            [maxstrain,Q,Chi,dChidz,phi0,dphi0dz,z] = NewtSolve(Pcal,Lcal,Scal,Tcal,Phi,zpts(i));
        end

         if isnan(dphi0dz)
            return
         end

      error =   norm( phi0-phi0prev,1)/norm(phi0,1);

        %Checks whether the solution has converged
        if  error <tol
            Converged=1;
            break;
        end
                
                % Interpolate current porosity profile to finer mesh for next iteration
                if i~=length(zpts)
                  phi0prev=interp1(z,phi0,linspace(0,1,zpts(i+1)),'pchip');
                end

            

    end
    
    %Displays warning if solution hasn't converged
    if Converged==0
        disp("Solution hasn't converged to within tol")
    end


end