function [Chi,dChidz,d2Chidz2] = MembraneProfApproxSol(Pcal,Lcal,Scal,Tcal,z)
%Function which returns the boundary layer approximation to the solution
%for the membrane equation

%Initalises varibales
G=1/3;
Mpl6G = (Scal+6*G);
alpha = (Mpl6G)/Scal;
beta= 1/Scal;
delta=Tcal^(1/4);
omega=(alpha/4)^(1/4)/delta;

%Boundary constants
A= -beta*(Pcal+Lcal)/alpha;
C =- beta*Pcal/alpha;
B=A+beta*Lcal/(omega*alpha);
D = C+beta*Lcal/(alpha*omega);

%Profile
Chi= exp(-omega*z).*(A*cos(omega*z)+B*sin(omega*z)) +exp(-omega*(1-z)).*(C*cos(omega*(1-z))+D*sin(omega*(1-z))) +beta*(Pcal+Lcal*(1-z))/alpha;

%First derivative of profile
dChidz= -omega*exp(-omega*z).*(A*cos(omega*z)+B*sin(omega*z)) + omega*exp(-omega*z).*(-A*sin(omega*z)+B*cos(omega*z)) ...
    +omega*exp(-omega*(1-z)).*(C*cos(omega*(1-z))+D*sin(omega*(1-z))) -omega*exp(-omega*(1-z)).*(-C*sin(omega*(1-z))+D*cos(omega*(1-z)))...
    -beta*Lcal/alpha;

%Second derivative of profile
d2Chidz2= omega^2*exp(-omega*z).*(A*cos(omega*z)+B*sin(omega*z)) -2*omega^2*exp(-omega*z).*(-A*sin(omega*z)+B*cos(omega*z)) ...
    + omega^2*exp(-omega*z).*(-A*cos(omega*z)-B*sin(omega*z)) ...
    +omega^2*exp(-omega*(1-z)).*(C*cos(omega*(1-z))+D*sin(omega*(1-z)))   -2*omega^2*exp(-omega*(1-z)).*(-C*sin(omega*(1-z))+D*cos(omega*(1-z)))...
    +omega^2*exp(-omega*(1-z)).*(-C*cos(omega*(1-z))-D*sin(omega*(1-z)));



end