function [data,BEAMLINE] = scanFunctionF2SingleBunchSim(simparams)

% add paths
addpath('simulationHelperFunctions')%for change linacphases func and other helper functions
% addpath(genpath('/Users/cemma/Documents/Work/FACET-II/Lucretia_sims/facet2-lattice'))% for loadfacetwakefields
addpath(genpath('/home/mstobbe/matlab_scripts/facet2-lattice-master'))% for loadfacetwakefields


% Load the latest electron beam for the 2nC single bunch case - Need to download this from the facet2 github repository
% load '/Users/cemma/Documents/Work/FACET-II/Lucretia_sims/LucretiaLatestGit/facet2-lattice/Lucretia/beams/FACET2e_op'
load '/home/mstobbe/matlab_scripts/facet2-lattice-master/Lucretia/beams/FACET2e_op.mat'

% load the latest lattice - this loads the entire 1 km lattice - Need to download this from the facet2 github repository
% load '/Users/cemma/Documents/Work/FACET-II/Lucretia_sims/LucretiaLatestGit/facet2-lattice/Lucretia/models/FACET2e/FACET2e'
load '/home/mstobbe/matlab_scripts/facet2-lattice-master/Lucretia/models/FACET2e/FACET2e.mat'

ISRon=0; % 1= Incoherent Sync. Rad. on, 0= ISR off
CSRon=0; % 1= Coherent Sync. Rad. on, 0= CSR off
LSCon=0; % 1= Longitudinal space charge on, 0= LSC off
% Decimate beam?
decBeam=100; % use integer>1 
             
% -- Tracking options
% Define indices
loadFACETLatticeIndices

% -- Add wakefields
BEAMLINE{194}.Wakes = [0 0];% This structure had wakes set to [4 0] which made the tracking crash
loadFacetWakefields
% Collective effects
%-- Sync. Radiation
eleSR=findcells(BEAMLINE,'Class','SBEN');
for ij=eleSR
  BEAMLINE{ij}.TrackFlag.SynRad=2*ISRon;
  BEAMLINE{ij}.TrackFlag.CSR=-1*CSRon;
  BEAMLINE{ij}.TrackFlag.CSR_SmoothFactor=0;
  BEAMLINE{ij}.TrackFlag.CSR_DriftSplit=25;
  BEAMLINE{ij}.TrackFlag.Split=25;
  % Turn off 2d CSR
        BEAMLINE{ij}.TrackFlag.CSR_2D = 0;
        BEAMLINE{ij}.TrackFlag.CSR_USEGPU = 0;
end

% -- Space charge
for ij=findcells(BEAMLINE,'TrackFlag')
  BEAMLINE{ij}.TrackFlag.LSC=LSCon;
  BEAMLINE{ij}.TrackFlag.LSC_storeData=0;
  % Set NBPM on LCAV elements to ensure 0.1m drift sections for
  % application of LSC
  if strcmp(BEAMLINE{ij}.Class,'LCAV')
    BEAMLINE{ij}.NBPM=LSCon*BEAMLINE{ij}.L/0.1;
    BEAMLINE{ij}.GetSBPMData=LSCon;
    BEAMLINE{ij}.GetInstData=LSCon;
  end
end

% Set initial beam parameters
Initial.Q=Initial.Q*(1+simparams.qi); % Initial charge offset
Beam.Bunch.x(1,:)=Beam.Bunch.x(1,:)+simparams.dx;% Initial offset in x
Beam.Bunch.x(3,:)=Beam.Bunch.x(3,:)+simparams.dy;% Initial offset in y

% --- Set Linac Phases / degrees off-crest and amplitudes
[BEAMLINE]=changeLinacPhases(BEAMLINE,simparams);

% Decimate beam
if decBeam>1
  decBeam=floor(decBeam);
  Beam.Bunch.x=Beam.Bunch.x(:,1:decBeam:end);
  Beam.Bunch.stop=Beam.Bunch.stop(1:decBeam:end);
  Beam.Bunch.Q=Beam.Bunch.Q(1:decBeam:end).*decBeam;
  disp('Beam decimated')
end   
plotSpotSizeFromAtoB(BEAMLINE,Initial,profMontcavS10,istart);

%% Track from the end of the injector sim (before the LH chicane) to the LH center
[~,beam]=TrackThru(istart,iLH,Beam,1,1,0); 
% Add Sinusoidal modulation to the beam
[beam,~,~,~]=addLHmodulation(beam,simparams.Elaser,0);
% Save the beam at LH center after applying LH modulation
data(1).beam = beam;
%% Continue tracking form the Laser Heater to the end of BC20
% Track to end of each bending section and re-center beam in each case to
% take care of phase-slip with respect to RF and orbit excursions due to SR
% energy losses. Also remove linear dispersion induced by CSR effects.
% In reality this is done with RF phasing, orbit feedbacks and beam tuning.
% Track through the laser heater

% FFStart = 1517 name = 'BEGFF20'
% FFend = 1568 name = 'ENDFF20'
ele=[iLH BC11END BC14END BC20END];
for ij=1:length(ele)-1
  [~,beam]=TrackThru(ele(ij),ele(ij+1),beam,1,1);
  disp(BEAMLINE{ele(ij)}.Name)
  beam = reCenterBeamFromRFPhaseSlip(beam);
  data(ij+1).beam = beam;
  beamImage(beam);
  [nx,ny,~] = GetNEmit90FromBeam( beam, 1);
  disp(sprintf(['Emit at ',BEAMLINE{ele(ij+1)}.Name,' = ',num2str(nx),', ', num2str(ny)]))
end
%% Track from BC20 end to the IP and plot the output LPS at the IP
[~,data(length(ele)+1).beam]=TrackThru(ele(end),PENT,beam,1,1);

        ebins = [9.5:0.01:10.2];% GeV % The resolution of ~0.5 MeV/pix comes from the pixel size/dispersion at TCAV
        zbins = ([-150:0.5:100]*1e-6);% m % The resolution of ~0.25 um/pix comes from the pixel size/streak in um/um
        %[~,~,~,data.I,data.Eprof]=MakeBeamLPS(beam,zbins,ebins,1);        
        [~,~,~,data(length(ele)+1).I,data(length(ele)+1).Eprof]=MakeBeamLPS(beam,zbins,ebins,1);        
end

