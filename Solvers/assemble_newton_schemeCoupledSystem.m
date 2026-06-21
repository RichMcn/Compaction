function [q, J] = assemble_newton_schemeCoupledSystem(v, z, Phi, Pcal, Lcal, Tcal, Scal, vphi)
%Assembles the Newton scheme for the full coupled system, without using approximations for Chi. In
%the paper this is used for benchmarking the solutions found via the
%approximations for Chi.

%Initialises variables
M = length(z);
dz = z(2) - z(1);
N = 2*M - 6;
q = zeros(N,1);
J = sparse(N,N);

%Helper functions
A = @(vi) Phi*vi/((1-vi)*(1-Phi));
Ap = @(vi) Phi/((1-vi)^2*(1-Phi));
App = @(vi) 2*Phi/((1-vi).^3*(1-Phi));
B = @(vi) 3*vi^3/(1-vi);
Bp = @(vi) 3*vi^2*(3-2*vi)/(1-vi)^2;
Bpp = @(vi) 6*vi*(3-3*vi+vi^2)/(1-vi)^3;

%-------------------------- Formation of the residual ------------------

%%— residual for φ-equations 
q(1) = A(v(1))*(1+v(M-1)/4)*(v(2)-2*v(1)+vphi)/dz^2 + Ap(v(1))*(1+v(M-1)/4)*( (v(2)-vphi)/(2*dz) )^2 + ...
    (4*A(v(1))-Bp(v(1))*(1+v(M-1)/4))*v(M-1)*(v(2)-vphi)/(4*dz^2) -B(v(1))*(1+v(M-1)/4)*v(M-1)/(2*dz^2) -4*B(v(1))*(v(M-1)/(2*dz))^2;

q(2) = A(v(2))*(1+v(M-1))*(v(3)-2*v(2)+v(1))/dz^2 +Ap(v(2))*(1+v(M-1))*( (v(3)-v(1))/(2*dz) )^2 + ...
    (4*A(v(2))-Bp(v(2))*(1+v(M-1)))*(v(M)-v(M-1)/4)*(v(3)-v(1))/(4*dz^2) - B(v(2))*(1+v(M-1))*(v(M)-7*v(M-1)/4)/dz^2 -4*B(v(2))*( (v(M)-v(M-1)/4)/(2*dz) )^2;

% qi, 2<i<M-3  
for i = 3 : M-4

    q(i) = A(v(i))*(1+v(i+M-3))*(v(i+1)-2*v(i)+v(i-1))/dz^2 + Ap(v(i))*(1+v(i+M-3))*( (v(i+1)-v(i-1))/(2*dz) )^2 + ...
        (4*A(v(i))-Bp(v(i))*(1+v(i+M-3)))*(v(i+M-2)-v(i+M-4))*(v(i+1)-v(i-1))/(4*dz^2) -B(v(i))*(1+v(i+M-3))*(v(i+M-2)-2*v(i+M-3)+v(i+M-4))/dz^2-4*B(v(i))*( (v(i+M-2)-v(i+M-4))/(2*dz) )^2;

end

q(M-3) = A(v(M-3))*(1+v(2*M-6))*(v(M-2)-2*v(M-3)+v(M-4))/dz^2 +Ap(v(M-3))*(1+v(2*M-6))*( (v(M-2)-v(M-4))/(2*dz) )^2 + ...
    (4*A(v(M-3))-Bp(v(M-3))*(1+v(2*M-6)))*(v(2*M-6)/4-v(2*M-7))*(v(M-2)-v(M-4))/(4*dz^2) - B(v(M-3))*(1+v(2*M-6))*(v(2*M-7)-7*v(2*M-6)/4)/dz^2 -4*B(v(M-3))*( (v(2*M-6)/4-v(2*M-7))/(2*dz) )^2;


q(M-2) = A(v(M-2))*(1+v(2*M-6)/4)*(Phi-2*v(M-2)+v(M-3))/dz^2 + Ap(v(M-2))*(1+v(2*M-6)/4)*( (Phi-v(M-3))/(2*dz) )^2 + ...
    (4*A(v(M-2))-Bp(v(M-2))*(1+v(2*M-6)/4))*(-v(2*M-6))*(Phi-v(M-3))/(4*dz^2) -B(v(M-2))*(1+v(2*M-6)/4)*v(2*M-6)/(2*dz^2) -4*B(v(M-2))*(-v(2*M-6)/(2*dz))^2;


q(M-1) = (Tcal/dz^4)*Scal*(v(M+1)-4*v(M)+5*v(M-1)) +(Scal+2)*v(M-1) - Pcal -Lcal*(1-z(3)) -2*((v(2)-Phi)/(1-Phi))/3;

q(M) =   (Tcal/dz^4)*Scal*(v(M+2) - 4*v(M+1) + 6*v(M) - 15*v(M-1)/4) +(Scal+2)*v(M)  - Pcal - Lcal*(1 - z(4))  - 2*( (v(3)-Phi)/(1-Phi)  )/3;

for i = M+1 : 2*M-8
  q(i) =      (Tcal/dz^4)*Scal*(v(i+2)-4*v(i+1)+6*v(i)-4*v(i-1)+ v(i-2)) +(Scal+2)*v(i) ...
    - Pcal - Lcal*(1 - z(i-(M-4)))    - 2*( (v(i-(M-3)) - Phi)/(1 - Phi) )/3;
end

q(2*M-7) =   (Tcal/dz^4)*Scal*( -15*v(2*M-6)/4 + 6*v(2*M-7) - 4*v(2*M-8) + v(2*M-9) ) ...
  +(Scal+2)*v(2*M-7)   - Pcal - Lcal*(1 - z(M-3))  - 2*( (v(M-4)-Phi)/(1-Phi)  )/3;


q(2*M-6) =   (Tcal/dz^4)*Scal*( v(2*M-8) - 4*v(2*M-7) + 5*v(2*M-6) ) ...
  +(Scal+2)*v(2*M-6) ...
  - Pcal - Lcal*(1 - z(M-2))   - 2*( (v(M-3)-Phi)/(1-Phi)  )/3;



% -------------------  Formation of the Jacobian --------------------

% — First row eqs 

J(1,1)   = -A(v(1))*(1+v(M-1)/4)*2/dz^2 +Ap(v(1))*(1+v(M-1)/4)*(v(2)-2*v(1)+vphi)/dz^2 ...
    +App(v(1))*(1+v(M-1)/4)*( (v(2)-vphi)/(2*dz) )^2 + ...
    (4*Ap(v(1))-Bpp(v(1))*(1+v(M-1)/4))*v(M-1)*(v(2)-vphi)/(4*dz^2) ...
    - Bp(v(1))*(1+v(M-1)/4)*v(M-1)/(2*dz^2) -4*Bp(v(1))*(v(M-1)/(2*dz))^2;

J(1,2)   = A(v(1))*(1+v(M-1)/4)/dz^2 + Ap(v(1))*(1+v(M-1)/4)*(v(2)-vphi)/(2*dz^2) ...
    +(4*A(v(1))-Bp(v(1))*(1+v(M-1)/4))*v(M-1)/(4*dz^2);


J(1,M-1) =  A(v(1))*(v(2)-2*v(1)+vphi)/(4*dz^2)+ Ap(v(1))*( (v(2)-vphi)/(2*dz) )^2/4+...
    (-Bp(v(1))/4)*v(M-1)*(v(2)-vphi)/(4*dz^2)+ (4*A(v(1))-Bp(v(1))*(1+v(M-1)/4))*(v(2)-vphi)/(4*dz^2)-...
    B(v(1))*v(M-1)/(8*dz^2)-B(v(1))*(1+v(M-1)/4)/(2*dz^2)-4*B(v(1))*v(M-1)/(2*dz^2);

% — Second row eqs 
J(2,1) = A(v(2))*(1+v(M-1))/dz^2  - Ap(v(2))*(1+v(M-1))*(v(3)-v(1))/(2*dz^2)-...
    (4*A(v(2))-Bp(v(2))*(1+v(M-1)))*v(M-1)/(4*dz^2);

J(2,2) = -2*A(v(2))*(1+v(M-1))/dz^2 + Ap(v(2))*(1+v(M-1))*(v(3)-2*v(2)+v(1))/dz^2 + App(v(2))*(1+v(M-1))*( (v(3)-v(1))/(2*dz) )^2 ...
   +(4*Ap(v(2))-Bpp(v(2))*(1+v(M-1)))*(v(3)-v(1))*v(M-1)/(4*dz^2) - Bp(v(2))*(1+v(M-1))*(v(M)-7*v(M-1)/4)/dz^2 - 4*Bp(v(2))*((v(M)-v(M-1)/4)/(2*dz))^2;


J(2,3) = A(v(2))*(1+v(M-1))/dz^2  + Ap(v(2))*(1+v(M-1))*(v(3)-v(1))/(2*dz^2)+...
    (4*A(v(2))-Bp(v(2))*(1+v(M-1)))*v(M-1)/(4*dz^2);


J(2,M-1) = A(v(2))*(v(3)-2*v(2)+v(1))/dz^2+Ap(v(2))*((v(3)-v(1))/(2*dz))^2 ...
    -(4*A(v(2))-Bp(v(2))*(1+v(M-1)))*(v(3)-v(1))/(16*dz^2)-...
    Bp(v(2))*(v(3)-v(1))*(v(M)-v(M-1)/4)/(4*dz^2) ...
    -B(v(2))*(1+v(M-1))*(-7/4)/dz^2-B(v(2))*(v(M)-7*v(M-1)/4)/dz^2+B(v(2))*(v(M)-v(M-1)/4)/(2*dz^2);

J(2,M) = (4*A(v(2))-Bp(v(2))*(1+v(M-1)))*(v(3)-v(1))/(4*dz^2)-B(v(2))*(1+v(M-1))/dz^2-4*B(v(2))*(v(M)-v(M-1)/4)/(2*dz^2);

% — interior φ-block eqs 
for i = 3 : M-4

  J(i,i-1)     = A(v(i))*(1+v(i+M-3))/dz^2 - Ap(v(i))*(1+v(i+M-3))*(v(i+1)-v(i-1))/(2*dz^2)...
      -(4*A(v(i))-Bp(v(i))*(1+v(i+M-3)))*(v(i+M-2)-v(i+M-4))/(4*dz^2);

  J(i,i)       = -2*A(v(i))*(1+v(i+M-3))/dz^2 + Ap(v(i))*(1+v(i+M-3))*(v(i+1)-2*v(i)+v(i-1))/dz^2 ...
      +App(v(i))*(1+v(i+M-3))*( (v(i+1)-v(i-1))/(2*dz) )^2 + (4*Ap(v(i))-Bpp(v(i))*(1+v(i+M-3)))*(v(i+M-2)-v(i+M-4))*(v(i+1)-v(i-1))/(4*dz^2)...
      -Bp(v(i))*(1+v(i+M-3))*(v(i+M-2)-2*v(i+M-3)+v(i+M-4))/dz^2 -4*Bp(v(i))*((v(i+M-2)-v(i+M-4))/(2*dz))^2;

  J(i,i+1)     = A(v(i))*(1+v(i+M-3))/dz^2 + Ap(v(i))*(1+v(i+M-3))*(v(i+1)-v(i-1))/(2*dz^2)...
      +(4*A(v(i))-Bp(v(i))*(1+v(i+M-3)))*(v(i+M-2)-v(i+M-4))/(4*dz^2);

  J(i,i+(M-4)) = -(4*A(v(i))-Bp(v(i))*(1+v(i+M-3)))*(v(i+1)-v(i-1))/(4*dz^2) ...
      -B(v(i))*(1+v(i+M-3))/dz^2 +4*B(v(i))*((v(i+M-2)-v(i+M-4))/(2*dz^2));

  J(i,i+(M-3)) =  A(v(i))*(v(i+1)-2*v(i)+v(i-1))/(dz^2) +Ap(v(i))*((v(i+1)-v(i-1))/(2*dz))^2 ...
  -Bp(v(i))*(v(i+M-2)-v(i+M-4))*(v(i+1)-v(i-1))/(4*dz^2) +2*B(v(i))*(1+v(i+M-3))/dz^2 ...
  -B(v(i))*(v(i+M-2)-2*v(i+M-3)+v(i+M-4))/dz^2;

  J(i,i+(M-2)) = (4*A(v(i))-Bp(v(i))*(1+v(i+M-3)))*(v(i+1)-v(i-1))/(4*dz^2) ...
      -B(v(i))*(1+v(i+M-3))/dz^2 -4*B(v(i))*((v(i+M-2)-v(i+M-4))/(2*dz^2));

end

J(M-3,M-4) = A(v(M-3))*(1+v(2*M-6))/dz^2 -Ap(v(M-3))*(1+v(2*M-6))*(v(M-2)-v(M-4))/(2*dz^2) ...
    -(4*A(M-3)-Bp(v(M-3))*(1+v(2*M-6)))*(v(2*M-6)/4-v(2*M-7))/(4*dz^2);


J(M-3,M-3) = -2*A(v(M-3))*(1+v(2*M-6))/dz^2 + Ap(v(M-3))*(1+v(2*M-6))*(v(M-2)-2*v(M-3)+v(M-4))/dz^2 ...
+App(v(M-3))*(1+v(2*M-6))*( (v(M-2)-v(M-4))/(2*dz) )^2 ...
+(4*Ap(M-3)-Bpp(v(M-3))*(1+v(2*M-6)))*(v(2*M-6)/4-v(2*M-7))*(v(M-2)-v(M-4))/(4*dz^2)   ...
-Bp(v(M-3))*(1+v(2*M-6))*(-7*v(2*M-6)/4+v(2*M-7))/dz^2 - 4*Bp(v(M-3))*((v(2*M-6)/4-v(2*M-7))/(2*dz))^2;

J(M-3,M-2) = A(v(M-3))*(1+v(2*M-6))/dz^2 +Ap(v(M-3))*(1+v(2*M-6))*(v(M-2)-v(M-4))/(2*dz^2) ...
    +(4*A(M-3)-Bp(v(M-3))*(1+v(2*M-6)))*(v(2*M-6)/4-v(2*M-7))/(4*dz^2);


J(M-3,2*M-7) = -(4*A(v(M-3))-Bp(v(M-3))*(1+v(2*M-6)))*(v(M-2)-v(M-4))/(4*dz^2) -B(v(M-3))*(1+v(2*M-6))/dz^2 ...
    +4*B(v(M-3))*(v(2*M-6)/4-v(2*M-7))/(2*dz^2);


J(M-3,2*M-6) = A(v(M-3))*(v(M-2)-2*v(M-3)+v(M-4))/dz^2 +Ap(v(M-3))*((v(M-2)-v(M-4))/(2*dz))^2 ...
    +(4*A(v(M-3))-Bp(v(M-3))*(1+v(2*M-6)))*(v(M-2)-v(M-4))/(16*dz^2) ...
    - Bp(v(M-3))*(v(2*M-6)/4-v(2*M-7))*(v(M-2)-v(M-4))/(4*dz^2) ...
    -B(v(M-3))*(-7*v(2*M-6)/4+v(2*M-7))/dz^2 ...
    - B(v(M-3))*(1+v(2*M-6))*(-7/4)/dz^2-B(v(M-3))*(v(2*M-6)/4-v(2*M-7))/(2*dz^2);



% Row M-2 (i = M-2)
J(M-2,M-3) = A(v(M-2))*(1+v(2*M-6)/4)/dz^2 - Ap(v(M-2))*(1+v(2*M-6)/4)*(Phi-v(M-3))/(2*dz^2) ...
    + (4*A(v(M-2))-Bp(v(M-2))*(1+v(2*M-6)/4))*v(2*M-6)/(4*dz^2);


J(M-2,M-2) = -2*A(v(M-2))*(1+v(2*M-6)/4)/dz^2 +Ap(v(M-2))*(1+v(2*M-6)/4)*(Phi-2*v(M-2)+v(M-3))/dz^2 ...
    + App(v(M-2))*(1+v(2*M-6)/4)*( (Phi-v(M-3))/(2*dz) )^2  ...
    -(4*Ap(v(M-2))-Bpp(v(M-2))*(1+v(2*M-6)/4))*v(2*M-6)*(Phi-v(M-3))/(4*dz^2) ...
    - Bp(v(M-2))*(1+v(2*M-6)/4)*v(2*M-6)/(2*dz^2) -4*Bp(v(M-2))*(v(2*M-6)/(2*dz))^2;

J(M-2,2*M-6) = A(v(M-2))*(Phi-2*v(M-2)+v(M-3))/(4*dz^2) +Ap(v(M-2))*((Phi-v(M-3))/(2*dz))^2/4 ...
   - (-Bp(v(M-2))/4)*(v(2*M-6))*(Phi-v(M-3))/(4*dz^2) ...
   - (4*A(v(M-2))-Bp(v(M-2))*(1+v(2*M-6)/4))*(Phi-v(M-3))/(4*dz^2) ...
   -B(v(M-2))*v(2*M-6)/(8*dz^2) - B(v(M-2))*(1+v(2*M-6)/4)/(2*dz^2)-4*B(v(M-2))*v(2*M-6)/(2*dz^2);

% — F-block Jacobian eqns 

% Row M-1
J(M-1,2)     = -2/(3*(1-Phi));
J(M-1,M-1)   =  Tcal*Scal*5/dz^4+Scal+2;
J(M-1,M)     = -4*Tcal*Scal/dz^4;
J(M-1,M+1)   =  Tcal*Scal/dz^4;

% Row M
J(M,3)       = -2/(3*(1-Phi));
J(M,M-1)     = -15*Tcal*Scal/(4*dz^4);
J(M,M)       =  Tcal*Scal*6/dz^4+Scal+2;
J(M,M+1)     =  -4*Tcal*Scal/dz^4;
J(M,M+2)     =  Tcal*Scal/dz^4;

% Interior F rows
for i = M+1 : 2*M-8
  J(i,i-(M-3)) = -2/(3*(1-Phi));
  J(i,i-2)     =  Tcal*Scal/dz^4;
  J(i,i-1)     =  -4*Tcal*Scal/dz^4;
  J(i,i)       =  6*Tcal*Scal/dz^4+Scal+2;
  J(i,i+1)     =  -4*Tcal*Scal/dz^4;
  J(i,i+2)     =   Tcal*Scal/dz^4;
end

% Row 2M-7
J(2*M-7,M-4)   = -2/(3*(1-Phi));
J(2*M-7,2*M-9) =  Tcal*Scal/dz^4;
J(2*M-7,2*M-8) = -4*Tcal*Scal/dz^4;
J(2*M-7,2*M-7) =  6*Tcal*Scal/dz^4+Scal+2;
J(2*M-7,2*M-6) = -15*Tcal*Scal/(4*dz^4);

% Row 2M-6
J(2*M-6,M-3)   = -2/(3*(1-Phi));
J(2*M-6,2*M-8) =  Tcal*Scal/dz^4;
J(2*M-6,2*M-7) =  -4*Tcal*Scal/dz^4;
J(2*M-6,2*M-6) =  5*Tcal*Scal/dz^4+Scal+2;

end
