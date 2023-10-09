function [beamOutAbove, beamOutBelow] = collimate(beam, E, width)
    
    beamOutAbove = struct();
    beamOutAbove.BunchInterval = beam.BunchInterval;
    beamOutAbove.Bunch.x = [];
    beamOutAbove.Bunch.Q = [];
    beamOutAbove.Bunch.stop = [];

    beamOutBelow = struct();
    beamOutBelow.BunchInterval = beam.BunchInterval;
    beamOutBelow.Bunch.x = [];
    beamOutBelow.Bunch.Q = [];
    beamOutBelow.Bunch.stop = [];

    energies = beam.Bunch.x(6,:);

    for ii = 1:length(energies)
        if (energies(ii) > E+width)
            beamOutAbove.Bunch.x = [beamOutAbove.Bunch.x beam.Bunch.x(:,ii)];
            beamOutAbove.Bunch.Q = [beamOutAbove.Bunch.Q beam.Bunch.Q(ii)];
            beamOutAbove.Bunch.stop = [beamOutAbove.Bunch.stop beam.Bunch.stop(ii)];

        elseif(energies(ii) < E-width)
            beamOutBelow.Bunch.x = [beamOutBelow.Bunch.x beam.Bunch.x(:,ii)];
            beamOutBelow.Bunch.Q = [beamOutBelow.Bunch.Q beam.Bunch.Q(ii)];
            beamOutBelow.Bunch.stop = [beamOutBelow.Bunch.stop beam.Bunch.stop(ii)];

        end
    end

end