function [pressurecurve]=LineFinderDarcyBackwards(Scal,Tcal,Pcal,Lcal,Phi)
% LineFinderDarcyBackwardsNewt
% Determines the curve in parameter space where the flow behavior
% transitions from Darcy-like (undeformed porous medium) to a deformation 
% influenced flow rate. 
%
% This is done by fixing one flow parameter (P or L),
% and varying the other. For each stiffness ratio value, the flow rate is
% computed across a range of pressures or gravitational parameters until the
% deviation from Darcy flow exceeds a tolerance.
%
% Inputs:
%   Scal   - Vector of stiffness ratios
%   Tcal   - Membrane bending
%   Pcal   - Vector or scalar of pressure values
%   Lcal   - Vector or scalar of gravity parameter values
%   Phi  - Initial porosity
%
% Output:
%   pressurecurve - Vector of critical values in P or L space at which flow deviates
%                   from Darcy prediction, for each value of the varied material parameter.


addpath('../Solvers')

%Tolerance for calculator
tol=10^(-3);

% Ensure exactly one Pcal or Lcal should be a scalar value, with the other a vector.
l_Ds=length(Scal);
l_Ps=length(Pcal);
l_Ls=length(Lcal);
if (l_Ds  <2 )
    error('Scal  should be a vector.')
end
if (l_Ps > 1 && l_Ls > 1)
    error('Only one of Ps or Ls should be a vector.')
end
dropfurther=0;



    T=Tcal;

%Initialises vector of indices for where the line is located
Pindices=nan(1,length(Scal));

if l_Ps>1
    pressures=Pcal;
    L=Lcal;
    darcpress= pressures+L;
else
    pressures=Lcal;
    P=Pcal;
    darcpress = pressures+P;
end

% Compute Darcy flow rate for all pressures once (used as reference)
Q_darcy = -Phi^3 * pi * darcpress / (1 - Phi)^2;

%The code calculated the location of the line backwards in parameter space
istart =length(Scal);
jstart=5;
finish=0;

   for i=istart:-1:1

       if (jstart==(length(pressures)))||finish
           break;
       end

       if jstart<1
           jstart=1;
       end

       j=jstart;

       while j<length(pressures)

           if j<1
               break;
           end

           if l_Ds>1
               D=Scal(i);
           else

              

           end

           %Updated pressure
           if l_Ps>1
               P=Pcal(j);
           else
               L=Lcal(j);
           end

        %Calculates flow rate
        [maxstrain,Q,~,~,~,~,~] = SolveUntilTol(P,L,D,T,Phi,1e-4,1);

        if isnan(Q)
            break
        end

        if maxstrain>0.25
            finish=1;
            break;
        end

        % Relative deviation from Darcy flow
        rel_error = abs((Q - Q_darcy(j)) / Q);

        % Flow has deviated from Darcy prediction beyond tolerance.
        % Record the current material parameter index and the corresponding pressure
        % value just *before* this deviation.
        if   rel_error >tol

            if j<2
                %If the value pressures(1) for Scal(i) is already
                %deviating from the Darcy velocity, we break the inner loop
                %and move to the next value in Scal. 
                break;
            end
            
            %Having immediately started on another material parameter we
            %may need to drop down further in order to find the location of
            %the line
            if dropfurther==1
                jstart=j-2;
                j=j-2;
                continue;
            end

            %Records the location of the line, and initialises the index of
            %the first test location for the pressure in the next bunch of
            %iterations to find the line. Set a variable dropfurther to 1
            %which will make the code drop further if we are still
            %immediately above the line
            Pindices(i) = j-1;
            jstart=j-2;
            dropfurther=1;
                
            break;

        end
        
        %If we are not above the line we get to this part of the code which
        %allows us to go up and find the line
        dropfurther=0;

        %Prints to the screen where we are in the space so that we can
        %understand how much more time the code might run 
        [i,j]

        j=j+1;

       end

   end

   %Enters the curve into the vector pressurecurve
   pressurecurve = nan(size(Scal));
   valid = ~isnan(Pindices);
   pressurecurve(valid) = pressures(Pindices(valid));

end


