function [mxstr,Qz0] = AsympFR(Phi,varphi,Pcal,Scal, Tcal)
%Function which returns the asymptotic solution for the flow rate

B = varphi/Phi;
gamma = 2*Pcal/(Scal+2);
Omega = ((Scal+2)/(4*Scal*Tcal))^(1/4);

z = linspace(0,1,1000);

[Chi,~,~] = MembraneProfApproxSol(Pcal,0,Scal,Tcal,z);

[phi0] = Asympphi(Pcal,Scal,Phi,varphi,z,Omega);

mxstr = (phi0-Phi)/(1-Phi) - 2*Chi;

mxstr = max(abs(mxstr));

 Qz0 = pi*Phi^3*B*(1+B)*Pcal*(1+3*gamma*(2+B^2/(1+B))-Phi*(2*(1-B^3)/(3*(1-B^2))-1 ...
     +(1-Phi)*B*log(B^2)/(2*(1-B)) )/(1-Phi) )  /(2*(1-Phi))   ;  



end