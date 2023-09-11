% This script allows you to run a Lucretia simulation tracking the beam 
% from the injector to the FACET-II IP.
% C Emma, Jul 2023

% User needs to set the simulation input parameters (see simparams below) 
% and the path to the FACET-II lattice and FACET-II e-beam at the start point. 
% The paths for the beam and lattice needs to be set in
% scanFunctionF2SingleBunchSim.m
% The lattice and e-beam files can both be found on github
% https://github.com/slaclab/facet2-lattice/Lucretia.

% Simulation input parameters


simparams = struct('P1',[+1 0.0], 'P2',[-1.65 0], 'V1',[5e-2 0.0], 'V2',[3e-2 0.0], ...
              'qi',[0 0.0], 'dx',[0 0], 'dy',[0 0],'Elaser',[0e-3 0]); % Center vals, deltas

% Notes on L1 and L2 phase values
% L1 phase at 0 delta = -20.5 deg
% L2 phase at 0 delta = -38.25 deg, note full comp at delta = -1.65 deg, moderate comp at delta = -1.8 deg

N = 1; % Number of simulations to run (if you want to do a scan with random sampling of input parameters)
rng('shuffle');

for ij = 1:N
    tic    
        vals = structfun(@samplevals,simparams,'UniformOutput',false); % Randomly samples parameter space
        data = scanFunctionF2SingleBunchSim(vals);%Without E fb    
        data(end).scanvals = vals;
        save('simulationData','data','-v7'); % Saves the simulation data
    toc
end
%% Example of looking at the beam simulation data after a run
simulationData = importdata('simulationData.mat');
beamImage(simulationData(end).beam)
%% Function to randomly sample number between two values
function val = samplevals(in)
low = in(1) - in(2);
high = in(1) + in(2);
val = (high-low)*rand()+low;
end
