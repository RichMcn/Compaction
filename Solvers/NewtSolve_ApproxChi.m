function  [maxstrain,Q,Chi,dChidz,phi0,dphi0dz,z] = NewtSolve_ApproxChi(Pcal,Lcal,Scal,Tcal,Phi,zpts)
% Solves for phi_0(z) using a Newton-Raphson method using the approximated
% Chi(z) from the simplification of the coupled system and bounday layer
% solution for Chi.
% Inputs are:
% Pcal = pressure 
% Lcal = gravitational parameter
% Scal = stiffness ratio
% Tcal = bending parameter
% Phi = porosity of unstrained medium
% zpts = the number of evenly spaced gridpoints on the interval [0,1], on which the solution will
% be found

% Convergence tolerance for the Newton-Raphson iteration
tol=10^(-6);

% Compute boundary condition phi_0(0) at the outlet based on pressure and initial porosity
varphi = Phi/(1+9*(1-Phi)*(Pcal+Lcal)/4);


% Initial guess: use leading order asymptotic solution
% This typically falls within the Newton basin for small strain regimes
% Initialise mesh
z=linspace(0,1,zpts)';
phiGuess  = sqrt(varphi^2+(Phi^2-varphi^2)*z);
dz=z(2)-z(1);
var=[phiGuess(2:end-1)];   

% Initialize previous iterate for convergence check and residual vector
q=ones(zpts,1);

% Computes approximate Chi(z) and its derivatives.
% Aborts if Chi(z) is undefined, e.g.,is Tcal is too small
[Chi,dChidz,d2Fdz2] = MembraneProfApproxSol(Pcal,Lcal,Scal,Tcal,z);
if isnan(Chi)
    return
end

loopvar=1;
resi=1;

      while (resi>tol)||(norm(q,inf)>tol)
             % Assemble Newton step. For debugging, uncomment resid to print residual.
             [q, J] = assemble_newton_schemeApproxChi(var, z, Phi, Chi, dChidz,d2Fdz2, varphi);
             varprev=var;

             
                 var = var -  J\ q;

                 %remove colon to print
                resi = (norm((var - varprev), inf)/norm(var,inf));
                 
                 loopvar=loopvar+1;

                 if loopvar==100
                     break
                 end

      end

      if loopvar==100
          if resi>1
                phi0=nan;
                  Q=nan;
                    maxstrain=nan;
                  dphi0dz=nan;
              return
          else
                sprintf("Error with convergence")
                phi0=nan;
                  Q=nan;
                    maxstrain=nan;
                  dphi0dz=nan;
                  return
          end
      end

% Reconstruct full phi_0 and Chi profiles, and computes first derivatives 
% using finite differences + interpolation
phi0 = [varphi,var(1:zpts-2)',Phi];
dphi0dz1 = diff(phi0)'./diff(z);


if any(isnan(dphi0dz1))
    Q=nan;
    maxstrain=nan;
    dphi0dz=nan;
    return
end
dphi0dz = [(-3*varphi/2+2*phi0(2)-phi0(3)/2)/dz,(phi0(3:end)-phi0(1:end-2))/(2*dz),...
    -(-3*Phi/2+2*phi0(end-1)-phi0(end-2)/2)/dz ];


% Qvec should be constant (representing a uniform flow rate), but small numerical noise may be present.
% We compute Q as the mean of the central portion of Qvec to suppress boundary artifacts.
Qvec =-4*pi*Phi*phi0.*(1+Chi').^4.*dphi0dz./(9*(1-Phi)^2.*(1-phi0))+4*pi*phi0.^3.*(1+Chi').^4.*dChidz'./(3*(1-phi0).*(1-Phi));
Q = mean(Qvec(round(zpts/4):round(3*zpts/4) ));


%Calculated the maximum component \varepsilon_{zz} of strain
strains = (phi0-Phi)./(1-Phi) -2*Chi';
maxstrain = max(abs(strains));

end