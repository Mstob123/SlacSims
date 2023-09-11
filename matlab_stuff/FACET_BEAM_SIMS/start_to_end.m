clear all;

c0=299792458;
mE=9.1094e-31;
qE=1.6022e-19;
e0 = 8.8542e-12;

pDens = 1e16*1e6;

load '/Users/mwstobbe/Documents/SlacSims/matlab_stuff/facet2-lattice-master/Lucretia/models/FACET2e/FACET2e.mat'

addpath('FACET_BEAM_SIMS/simulationHelperFunctions')

oneBunch = load('/Users/mwstobbe/Documents/SlacSims/matlab_stuff/facet2-lattice-master/Lucretia/beams/FACET2e_op.mat');
twoBunch = load('/Users/mwstobbe/Documents/SlacSims/matlab_stuff/facet2-lattice-master/Lucretia/beams/FACET2e_2bunch_op.mat');

initialBeam = oneBunch.Beam;
initialParams = oneBunch.Initial;

twoBunchBeam = twoBunch.Beam;

simparams = struct('P1',[+1 0.0], 'P2',[-1.65 0], 'V1',[5e-2 0.0], 'V2',[3e-2 0.0], ...
              'qi',[0 0.0], 'dx',[0 0], 'dy',[0 0],'Elaser',[0e-3 0]);

vals = structfun(@samplevals,simparams,'UniformOutput',false);

FFstart = findcells(BEAMLINE,'Name','BEGFF20')
FFend = findcells(BEAMLINE,'Name','ENDFF20')
PENT = findcells(BEAMLINE,'Name','PENT')

beamImage(initialBeam)

tic
[data, BEAMLINE] = singleBunchTransport(vals,FFstart,initialBeam,BEAMLINE);
toc

beam = data(end).beam;
beamImage(beam)

[betax,betay] = findBeta(beam)

optimize = 1;

if (optimize == true)
    quadEle = findcells(BEAMLINE,'Class','QUAD', FFstart, FFend);

    gamavg = mean(beam.Bunch.x(6,:)*1000/0.511);
    betaMatch = c0/qE*sqrt(e0*mE)*sqrt(2*gamavg/pDens)

    wrappedBetaMin = @(v) betaMin(BEAMLINE,beam,quadEle,FFstart,PENT,betaMatch,v);

    options = optimset('Display','iter');
    for ii = 1:length(quadEle)/2
        x0(ii) = BEAMLINE{quadEle(2*ii)}.B;
    end

    x0

    x = fminsearch(wrappedBetaMin,x0, options)

    for ii = 1:2:length(quadEle)
        BEAMLINE{quadEle(ii)}.B = x((ii+1)/2);
        BEAMLINE{quadEle(ii+1)}.B = x((ii+1)/2);
    end

    [~,beamOut] = TrackThru(FFstart,PENT,beam,1,1);

else

    [~,beamOut] = TrackThru(FFstart,PENT,beam,1,1);

end

beamImage(beamOut)

fileOutd = '/Users/mwstobbe/Documents/SlacSims/matlab_stuff/driver_test.txt';
fileOutw = '/Users/mwstobbe/Documents/SlacSims/matlab_stuff/witness_test.txt';

Luc2FBPICtxt(beamOut,fileOutd)
Luc2FBPICtxt(beamOut,fileOutw)


%{
conda_environment = 'masontest'; % Replace with your Conda environment name
python_script = '/home/mstobbe/FBPIC_Sims/pwaSim.py'; % Replace with the path to your Python script
arguments = [fileOutd,' ',fileOutw] % Replace with your desired command line arguments

command = sprintf('conda run -pwd /home/mstobbe/FBPIC_Sims/ -n %s python3.11 %s %s', conda_environment, python_script, arguments);
system(command);
%}


%% Functions
function val = samplevals(in)
    low = in(1) - in(2);
    high = in(1) + in(2);
    val = (high-low)*rand()+low;
end


function match = betaMin(BEAMLINE, beamIn, quadEle, start, fin, betaMatch, v)

    for ii = 1:2:length(quadEle)
        BEAMLINE{quadEle(ii)}.B = v((ii+1)/2);
        BEAMLINE{quadEle(ii+1)}.B = v((ii+1)/2);
    end

    
    [~,beam] = TrackThru(start, fin, beamIn,1,1);

    [betax,betay] = findBeta(beam);

    match = sqrt((betax-betaMatch)^2+(betay-betaMatch)^2);

end

function [betax,betay] = findBeta(beam)

    x = beam.Bunch.x(1,:);
    xPrime = beam.Bunch.x(2,:);

    x_xPrime = x.*xPrime;

    x2_mean = mean(x.^2)-mean(x)^2;
    xPrime2_mean = mean(xPrime.^2)-mean(xPrime)^2;
    x_xPrime_mean2 = (mean(x_xPrime)-mean(x)*mean(xPrime))^2;

    emitx = sqrt(x2_mean*xPrime2_mean-x_xPrime_mean2)*mean(beam.Bunch.x(6,:)*1000/0.511);
    stdx = std(x);

    betax = stdx^2/emitx;

    y = beam.Bunch.x(3,:);
    yPrime = beam.Bunch.x(4,:);

    y_yPrime = y.*yPrime;

    y2_mean = mean(y.^2)-mean(y)^2;
    yPrime2_mean = mean(yPrime.^2)-mean(yPrime)^2;
    y_yPrime_mean2 = (mean(y_yPrime)-mean(y)*mean(yPrime))^2;

    emity = sqrt(y2_mean*yPrime2_mean-y_yPrime_mean2)*mean(beam.Bunch.x(6,:)*1000/0.511);
    stdy = std(y);

    betay = stdy^2/emity;

end



    






