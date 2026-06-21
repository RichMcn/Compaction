# Compaction
Code for reproducing figures from "Compaction in deformable porous media cylinder with elastic boundaries"

The code for creating all of the figures are in the main folder. Each one is a script which can be easily adapted.
The five model parameters are called Pcal, Lcal, Scal, Tcal, and Phi throughout.
PlotFig5... which calculates the curves for the regime map will take several hours to run.
In the main folder there is also a VelLoop script which can be used to sweep over various values of the parameters 
to find flow rates.
All of the functions used to solve the equations are in the Solvers script.
All of the functions which find the curves for the regime map are in the LineFinders folder.
