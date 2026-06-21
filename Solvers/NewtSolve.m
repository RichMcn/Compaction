function  [maxstrain,Q,Chi,dChidz,phi0,dphi0dz,z] = NewtSolve(Pcal,Lcal,Scal,Tcal,Phi,zpts)
% Solves the coupled nonlinear system for Chi(z) and phi_0(z) using a Newton-Raphson scheme
% Inputs are:
% Pcal = imposed fluid pressure 
% Lcal = gravitational pressure parameter
% Scal = ratio of stiffnesses
% Tcal = bending parameter
% Phi = porosity of unstrained medium
% zpts = the number of evenly spaced gridpoints on the interval [0,1], on which the solution will
% be found

% Convergence tolerance for the Newton-Raphson iteration
tol=10^(-6);

% Compute boundary condition phi_0(0) at the outlet based on pressure and initial porosity
varphi = Phi/(1+9*(1-Phi)*(Pcal+Lcal)/4);
z=linspace(0,1,zpts)';

%Initial guess for phi_0 and Chi(z). 
%phiGuess  = Phi*ones(zpts,1)/2;
phiGuess = sqrt(varphi^2+(Phi^2-varphi^2)*z);
FGuess = zeros(zpts,1);

var=[phiGuess(2:end-1);FGuess(3:end-2)];  

dz=z(2)-z(1);

% Scaling matrix to improve conditioning of the Jacobian system in certain parameter regimes
S = ones(size(var));
S(zpts-1:end) = dz^2;
SMatrix = spdiags(S,0,length(var),length(var));


resi=10*tol;

      while resi>tol
             % Assemble Newton step. For debugging, uncomment resid to print residual.
             % To slightly improve speed, use the unscaled update (var = var - J\q)
             [q, J] = assemble_newton_schemeCoupledSystem(var, z, Phi, Pcal, Lcal, Tcal, Scal, varphi);
             varprev=var;

              
            var = var - (SMatrix* J) \ (S .* q);

             resi = norm(var - varprev, inf)/norm(var, inf);
      end

    

% Reconstruct full phi_0 and Chi profiles, and computes first derivatives using finite differences + interpolation
phi0 = [varphi,var(1:zpts-2)',Phi];
Chi = [0,var(zpts-1)/4,var(zpts-1:2*zpts-6)',var(2*zpts-6)/4,0];
dFdz1 = diff(Chi)'./diff(z);
dChidz = [0,interp1(z(1:end-1)'+0.5*diff(z)',dFdz1,z(2:end-1),'pchip')',0 ];
dphi0dz1 = diff(phi0)'./diff(z);
dphi0dz = [(-3*varphi/2+2*phi0(2)-phi0(3)/2)/dz,interp1(z(1:end-1)'+0.5*diff(z)',dphi0dz1,z(2:end-1)','pchip'),-(-3*Phi/2+2*phi0(end-1)-phi0(end-2)/2)/dz ];

% Optional: plot porosity profile and deformed geometry (1+F vs. z)
% figure(1)
% plot(z,phi0)
% figure(2)
% plot(1+F,z)
% xlim([0.85,1.15])

% Qvec should be constant (representing a uniform flow rate), but small numerical noise may be present.
% We compute Q as the mean of the central portion of Qvec to suppress boundary artifacts.
Qvec =-4*pi*Phi*phi0.*(1+Chi).^4.*dphi0dz./(9*(1-Phi)^2.*(1-phi0))+4*pi*phi0.^3.*(1+Chi).^4.*dChidz./(3*(1-phi0).*(1-Phi));
Q = mean(Qvec(round(zpts/4):round(3*zpts/4) ));





%Calculated the maximum component \varepsilon_{zz} of strain
strains = -(phi0-Phi)./(1-phi0) -2*Chi;
maxstrain = max(abs(strains));


end