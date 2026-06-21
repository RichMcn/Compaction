function [phi0] = Asympphi(Pcal,Scal,Phi,varphi,z,Omega)
%Function which returns the asymptotic solution for the porosity

gamma = 2*Pcal/(Scal+2);

B = varphi/Phi;

phi0 =Phi* ( sqrt((1-B^2)*z+B^2) ...
    + gamma*(   ...
   3* ( B^3+(1-B^3)*z)./(2* sqrt((1-B^2)*z+B^2) )  - 3*exp(-Omega*z).*((1-B^2)*z+B^2).*(sin(Omega*z)+cos(Omega*z))/2 ...
    -3*exp(-Omega*(1-z)).*((1-B^2)*z+B^2).*(sin(Omega*(1-z))+cos(Omega*(1-z))/2 )    ...
    )    );

end