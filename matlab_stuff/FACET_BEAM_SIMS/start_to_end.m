function data = start_to_end(match_optimize, p_interact, plasma_density, L2_phase_offset)

close all;
data = struct();

c0=299792458;
mE=9.1094e-31;
qE=1.6022e-19;
e0 = 8.8542e-12;

%pDens = 1e16*1e6;
pDens = plasma_density;

load './matlab_stuff/facet2-lattice-master/Lucretia/models/FACET2e/FACET2e.mat'

%addpath('./matlab_stuff/FACET_BEAM_SIMS/simulationHelperFunctions')

oneBunch = load('./matlab_stuff/facet2-lattice-master/Lucretia/beams/FACET2e_op.mat');
twoBunch = load('./matlab_stuff/facet2-lattice-master/Lucretia/beams/FACET2e_2bunch_op.mat');

initialBeam = oneBunch.Beam;
initialParams = oneBunch.Initial;

data(1).beam = initialBeam;
twoBunchBeam = twoBunch.Beam;

%simparams = struct('P1',[+1 0.0], 'P2',[-1.65 0], 'V1',[5e-2 0.0], 'V2',[3e-2 0.0], ...
  %  'qi',[0 0.0], 'dx',[0 0], 'dy',[0 0],'Elaser',[0e-3 0]);

%vals = structfun(@samplevals,simparams,'UniformOutput',false);
%disp(vals)

vals = struct('P1',1, 'P2',-1.65+L2_phase_offset, 'V1',5e-2, 'V2',3e-2, ...
    'qi',0, 'dx',0 , 'dy',0 ,'Elaser',0e-3);



FFstart = findcells(BEAMLINE,'Name','BEGFF20')
FFend = findcells(BEAMLINE,'Name','ENDFF20')
PENT = findcells(BEAMLINE,'Name','PENT')
PEXT = findcells(BEAMLINE,'Name','PEXT')
DTOTR = findcells(BEAMLINE,'Name','DTOTR')

beamImage(initialBeam)

tic
[t_data, BEAMLINE] = singleBunchTransport(vals,FFstart,initialBeam,BEAMLINE);
toc

for ii = 1:length(t_data)
    data(ii+1).beam = t_data(ii).beam;
end

beam = t_data(end).beam;
beamImage(beam)

[betax,betay] = findBeta(beam);

if (match_optimize == true)
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

    [~,beam] = TrackThru(FFstart,PENT,beam,1,1);

else

    [~,beam] = TrackThru(FFstart,PENT,beam,1,1);

end

data(length(data)+1).beam = beam;
beamImage(beam);

fileOutd = './matlab_stuff/driver_test.txt';
fileOutw = './matlab_stuff/witness_test.txt';

Luc2FBPICtxt(beam,fileOutd)
Luc2FBPICtxt(beam,fileOutw)

if (p_interact == true)
    conda_environment = 'Facet_Simulations'; % Replace with your Conda environment name
    python_script = '.\python_stuff\pwaSim.py'; % Replace with the path to your Python script
    arguments = [fileOutd,' ',fileOutw];
    command = ['C:\Users\Mason\anaconda3\_conda run -n ', conda_environment, ' ', 'python '  python_script, ' ', arguments];

    [status, result] = system(command);
    disp(result)

    quadRange = linspace(-30,-36,10);

    folder = './diags/hdf5/';
    addpath(folder)

    files = dir(folder);
    files = files(~[files.isdir]);

    len = length(data);

    for ii = 1:numel(files)
        diag_beam = FBPIC_ReadElectron(folder, files(ii).name,'electrons_witness');
        lucBeam = YABP_FBPIC2YABP(diag_beam);
        data(ii+len).beam = lucBeam;
    end

    beam = data(end).beam;

    
else

    
    %beam.Bunch.x(:,5) = beamOut.Bunch.x(:,5) + pLength;

    %data(length(data)+1).beam = beam;

end

beam = trackThru_Spectrometer(beam,initialParams,PEXT,DTOTR,mean(beam.Bunch.x(6,:)));
data(length(data)+1).beam = beam;


%beamImage(beam)

save('data.mat','data')


end


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





    






