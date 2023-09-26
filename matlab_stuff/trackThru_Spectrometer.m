function beamOut = trackThru_Spectrometer(beamIn, Initial, iStart, iEnd, E)

% Tune the spectrometer, and track from the a point at the IP to the end of
% the spectrometer
% 
% Inputs
%  - beamIn  = the particle distribution in the standard Lucretia format
%  - Initial = the normal Lucretia initial structure (from the beginning of the BEAMLINE)
%  - iStart, iEnd = the starting and stoping element numbers
%                      - i.e. iStart = end of plasma, iEnd = DTOTR
%  - E = The energy that you want the quads to focus at (nominally the centroid energy of beamIn)
%
% Author - D. Storey, Sept 2023


% Must have the luctretia lattice loaded before running this function
global BEAMLINE PS

% % This will move all magnet strengths from the BEAMLINE struct to the PS structure
% if isempty(PS)
%     % Set up the lattice to use a PS
%     SetElementBlocks( 1, length(BEAMLINE) );
%     SetElementSlices( 1, length(BEAMLINE) );
%     SetIndependentPS( 1, length(BEAMLINE) );
% 
%     %Find all quads, and transfer magnet strengths to the PS
%     quads  = findcells(BEAMLINE,'Class','QUAD',1,length(BEAMLINE));
%     quadPS = [];
%     for iele=quads
%         quadPS(end+1) = BEAMLINE{iele}.PS;
%     end
%     quadPS = unique(quadPS);
%     MovePhysicsVarsToPS(quadPS)
% end

% Find the spectrometer quads in the lucretia deck
iQ0D = findcells(BEAMLINE,'Name','Q0D');
iQ1D = findcells(BEAMLINE,'Name','Q1D');
iQ2D = findcells(BEAMLINE,'Name','Q2D');

% Use the matching function to determine the quad settings to reimage
% between the element location of beamIn to the element of the diagnostic
% screen in the spectrometer
M12x = 0;    % this sets the reimaging condition
Ebend = 10;  % nominal dipole setting
tieQ02D = 1; % Normally these two are set to be equal

% Do the match
[Bdes] = match_SpectrometerSimple(Initial, E, iStart, iEnd, M12x, Ebend, tieQ02D);

% Set the quads to the strength you just found
PS(BEAMLINE{iQ0D(1)}.PS).Ampl = Bdes(1);
PS(BEAMLINE{iQ1D(1)}.PS).Ampl = Bdes(2);
PS(BEAMLINE{iQ2D(1)}.PS).Ampl = Bdes(3);


% Track the beam through
[~,beamOut]=TrackThru(iStart,iEnd,beamIn,1,1);

% View the output
beamImage(beamOut)

end