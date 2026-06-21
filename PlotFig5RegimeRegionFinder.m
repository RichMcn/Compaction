% PlotFig5RegimeRegionFinder
%
% This script computes and plots the various curves in the regime map
% figure. This script can take several hours to run depending on the chosen
% resolution
%
% Inputs / Parameter Setup:
%   - Phi: Porosity parameter (scalar)
%   - Scal: Vector of stiffness ratio
%   - Tcal: membrane bending parameter
%   - Pcal, Lcal: Pressure and gravity parameter vectors (one vector, one scalar)
%
%   - strainlimits: Array of strain thresholds to identify transition lines
%   - Plateautols: Array of tolerances defining proximity to plateau flow rate
%
% Outputs:
%   Figure showing the regime curves

addpath('LineFinders')

% Parameter sweep sizes
lengthpressuresweep = 12000;
lengthmaterialparamsweep = 12001;

% Porosity
Phi = 0.05;

% Define parameter spaces (only one of each pair should be vectors)
Pcal = logspace(-5, 6, lengthpressuresweep) ;   % Pressure
Lcal = 0;   % Gravity parameter 
Scal = logspace(-2, 7, lengthmaterialparamsweep);  % Flexural rigidity cleqr q
Tcal = 1e-07; % Membrane tension

% Strain limits for large strain detection curves
strainlimits = [ 0.01, 0.025, 0.05, 0.1,0.2];

% Tolerances for plateau flow rate detection
Plateautols = [2.5e-2, 1e-2, 1e-3];

% Preallocate matrices to store results
StrCurves = nan(length(strainlimits), max(length(Tcal), length(Scal)));
PlateauBottoms = nan(length(Plateautols), lengthmaterialparamsweep);
PlateauTops = nan(length(Plateautols), lengthpressuresweep);

% Clear temporary variables
clear lengthpressuresweep lengthmaterialparamsweep

% Compute Darcy transition curve
DarcyCurve = LineFinderDarcyBackwards(Scal, Tcal, Pcal, Lcal, Phi);

% Compute plateau flow rate curves for each tolerance level
for i = 1:length(Plateautols)
    % Bottom-up search for plateau transition line
    pressurecurve = LineFinderPlateauNewtBottom(Scal, Tcal, Pcal, Lcal, Phi, Plateautols(i));
    
    % Top-down vertical search for plateau transition line
    MaterialParamCurve = LineFinderPlateauNewtTopVert(Scal, Tcal, Pcal, Lcal, Phi, Plateautols(i));
    
    % Store results
   PlateauBottoms(i, :) = pressurecurve;
    PlateauTops(i, :) = MaterialParamCurve;
    
    % Clear temporary variables for memory management
    clear pressurecurve MaterialParamCurve
end

% Compute large strain exceedance curves for each strain limit
for i = 1:length(strainlimits)
    StrCurve = LineFinderLargeStrain(Scal, Tcal, Pcal, Lcal, Phi, strainlimits(i));
    StrCurves(i, :) = StrCurve;
    clear StrCurve
end

%Finds the curve, below which (below in terms of Scal) the contact condition fails
[CompStrainCurve] = LineFinderCompStr(Scal, Tcal, Pcal, Lcal, Phi);

%Finds the curve, northwest of which the porosity somewhere rises above the
%initial porosity
[compactioncurve]=LineFinderMaxphi(Scal,Tcal,Pcal,Lcal,Phi);

%If we are sweeping over Pcal, finds vertical curves which separate sub and
%super-Darcy. Also plots the solution
if ~isscalar(Pcal)
    %Finds vertical sections of the Darcy Curve
    [DarcyVertCurve] = LineFinderDarcyVert(Scal, Tcal, Pcal, Lcal, Phi,-10^(-3));
    [DarcyVertCurve2] = LineFinderDarcyVert(Scal, Tcal, Pcal, Lcal, Phi,10^(-3));
        PlotRegimeDiagramPcal
else
        PlotRegimeDiagramLcal
end

