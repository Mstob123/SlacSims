function [B] = match_Spectrometer(Initial, E, elezob, elezim, M12x, Ebend, tieQ02D)

% Do a match for various spectrometer settings, mainly the object plane
% 
% Inputs
%  - Initial - the normal Lucretia initial structure
%  - E - energy that you want to set the quads to
%  - zob, zim - the object and image planes
%  - M12x - the x plane M12. Assumes you always want to minimize M34
%  - Ebend - set nominal dispersion energy
%  - tieQ02D - =1 is fixing their values together, /=1 will allow them to differ
%       - note: you can't undo this without clearing the BEAMLINE and PS structures
%
% Author - D. Storey, June 2023
%             - simplified, Sept 2023


global BEAMLINE PS

% Reset the PS structure
for ips=1:length(PS)
  for iele=PS(ips).Element
    if GetTrueStrength(iele)==0
      BEAMLINE{iele}.B=0;
    else
      RenormalizePS(ips);
    end
    BEAMLINE{iele}.PS=0;
  end
end
PS = [];
SetIndependentPS( 1, length(BEAMLINE) );


if nargin<5
  error('Check input arguments');
end


% Get initial structure at object plane:
[~,T]=GetTwiss(1,elezob,Initial.x.Twiss,Initial.y.Twiss);
I=TwissToInitial(T,elezob,Initial);
I.Q=2e-9;
I.x.NEmit=3.0e-6; I.y.NEmit=3e-6;

% Find the spectrometer quads
quadele = findcells(BEAMLINE,'Class','QUAD',elezob, elezim);
quadPS = [];
for qe = quadele
    if ~ismember(BEAMLINE{qe}.PS, quadPS)
        quadPS = [quadPS BEAMLINE{qe}.PS];
    end
end
MovePhysicsVarsToPS( quadPS); % This just uses the power supply as the B value, rather then the beamline attribute
quadMAX = [239.7 386.8 233.2]/10;  % These are the real physical maxes of the power supplies

% Tie the first and last quad together
if exist('tieQ02D')
    if tieQ02D == 1
        PS(quadPS(1)).Element = [quadele([1 2 5 6])]; % Q0D and Q2D
        BEAMLINE{quadele(5)}.PS = quadPS(1);
        BEAMLINE{quadele(6)}.PS = quadPS(1);
        quadPS = quadPS(1:2);
        quadMAX = quadMAX(1:2);
    end
end


% Adjust element design energies
for ind = elezob:elezim 
    BEAMLINE{ind}.P = E;
end

if isempty('Ebend')
    bend=findcells(BEAMLINE,'Class','SBEN', elezob, elezim); % Find spectrometer magnet
    Bnom = 0.1001;
    for ii = 1:length(bend)
        BEAMLINE{bend(ii)}.B = Bnom*(Ebend/10);
    end
end


%% set up the match
M=Match;
M.initStruc=I; 
M.iInitial=elezob;
M.optim='fminsearch';

% Add spectrometer quads to match
for PSele = 1:length(quadPS) % the last three quads are the spectrometer qauds
    addVariable(M,'PS',quadPS(PSele),'Ampl',-quadMAX(PSele), quadMAX(PSele));
end

%Add constraints
M.addMatch(elezim,'R',M12x,1e-6,'12');
M.addMatch(elezim,'R',0,1e-6,'34');    % Assume you always want to minimize y-size

%Perform the matching
M.doMatch;
display(M)

[~,Twiss]=GetTwiss(elezob, elezim,I.x.Twiss,I.y.Twiss);
[~,R]=RmatAtoB(elezob, elezim)



B = [];
for qe = 1:2:length(quadele) 
    B(end+1) = GetTrueStrength(quadele(qe)) + GetTrueStrength(quadele(qe+1));
end
%TwissPlot(elezob, elezim, I, [1 0 1],.01);
