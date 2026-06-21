function [pressurecurve]=LineFinderMaxphi(Scal,Tcal,Pcal,Lcal,Phi)
% LineFinderDarcyBackwardsNewt
% Determines the boundary in parameter space where the flow behavior
% transitions from Darcy-like (undeformed porous medium) to deformation-dominated.
%
% This is done by fixing one flow parameter (Pcal or Lcal),
% and varying the other. For each stiffness ratio, the flow rate is
% computed across a range of pressures or gravitational parameters until the
% deviation from Darcy flow exceeds a tolerance.
%
% Inputs:
%   Scal   - Vector of stiffness ratios
%   Tcal   - Membrane bending parameter
%   Pcal   - Vector or scalar of pressure values
%   Lcal   - Vector or scalar of gravity parameter values
%   Phi  - Initial porosity
%
% Output:
%   pressurecurve - Vector of critical values in Pcal or Lcal space at which flow deviates
%                   from Darcy prediction, for each value of the varied material parameter.

addpath('../Solvers')


%Ps or Ls should be a scalar value, with the other a vector.
l_Ps=length(Pcal);
l_Ls=length(Lcal);
if (l_Ps > 1 && l_Ls > 1)
    error('Only one of Ps or Ls should be a vector.')
end
dropfurther=0;


%Initialises vector of indices for where the line is located
Pindices=nan(1,length(Scal));

if l_Ps>1
    pressures=Pcal;
    L=Lcal;
else
    pressures=Lcal;
    P=Pcal;
end


%The code calculated the location of the line backwards in parameter space
jstart=1;
finish =0;

   for i=1:length(Scal)

       if (jstart==(length(pressures)))
           break;
       end

       if jstart<1
           jstart=1;
       end

       if finish
           break;
       end

       j=jstart;

       while j<length(pressures)

           if j<1
               break;
           end

           if j==length(pressures)-1
               finish=1;
               break;
           end

               D=Scal(i);

           %Updated pressure
           if l_Ps>1
               P=Pcal(j);
           else
               L=Lcal(j);
           end

        %Calculates flow rate
        [~,~,~,~,phi0,~,~]  = SolveUntilTol(P,L,D,Tcal,Phi,1e-4,1);

        if isnan(phi0)
            break
        end

        % Relative deviation from Darcy flow
        maxphivar = (max(phi0)-Phi)/Phi;

        % Flow has deviated from Darcy prediction beyond tolerance.
        % Record the current material parameter index and the corresponding pressure
        % value just before this deviation.
        if   maxphivar >0

            if j<2
                %If the value pressures(1) for Scal(i) is already
                %deviating from the Darcy velocity, we break the inner loop
                %and move to the next value in sweep_param. 
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


