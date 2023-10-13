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


FFstart = findcells(BEAMLINE,'Name','BEGFF20');
FFend = findcells(BEAMLINE,'Name','ENDFF20');
eleFilg = findcells(BEAMLINE, 'Name', 'eleFILG');
PENT = findcells(BEAMLINE,'Name','PENT');
PEXT = findcells(BEAMLINE,'Name','PEXT');
DTOTR = findcells(BEAMLINE,'Name','DTOTR');

beamImage(initialBeam)

tic
[t_data, BEAMLINE] = singleBunchTransport(vals,FFstart,initialBeam,BEAMLINE);
toc

for ii = 1:length(t_data)
    data(ii+1).beam = t_data(ii).beam;
end

beam = t_data(end).beam;
f = figure('Name', 'Before Collimate Transverse');
f2 = figure('Name', 'Before Collimate Longitudinal');

beamImage(beam,[],[],[],[],[],[f, f2])

e_mean = mean(beam.Bunch.x(6,:))
[beam_drive,beam_wit] = collimate(beam, e_mean,0.003*e_mean);

f3 = figure('Name', 'Witness Transverse');
f4 = figure('Name', 'Witness Longitudinal');

f5 = figure('Name', 'Drive Transverse');
f6 = figure('Name', 'Drive Longitudinal');

beamImage(beam_wit,[],[],[],[],[],[f3, f4])
beamImage(beam_drive,[],[],[],[],[],[f5, f6])

%[betax,betay] = findBeta(beam_wit)

if (match_optimize == true)
    quadEle = findcells(BEAMLINE,'Class','QUAD', FFstart, FFend);

    gamavg = mean(beam_wit.Bunch.x(6,:)*1000/0.511);
    betaMatch = c0/qE*sqrt(e0*mE)*sqrt(2*gamavg/pDens);

    wrappedBetaMin = @(v) betaMin(BEAMLINE,beam_wit,quadEle,FFstart,PENT,betaMatch,v);

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

    [~,beam_wit] = TrackThru(FFstart,PENT,beam_wit,1,1);
    [~,beam_drive] = TrackThru(FFstart,PENT,beam_drive,1,1);

else

    [~,beam_wit] = TrackThru(FFstart,PENT,beam_wit,1,1);
    [~,beam_drive] = TrackThru(FFstart,PENT,beam_drive,1,1);

end

data(length(data)+1).beam = beam_wit;

fileOutd = './matlab_stuff/driver.txt';
fileOutw = './matlab_stuff/witness.txt';

Luc2FBPICtxt(beam_drive,fileOutd)
Luc2FBPICtxt(beam_wit,fileOutw)

ffff

if (p_interact == true)
    conda_environment = 'Facet_Simulations'; % Replace with your Conda environment name
    python_script = '.\python_stuff\pwaSim.py'; % Replace with the path to your Python script
    arguments = [fileOutd,' ',fileOutw, ' ', num2str(pDens)];
    command = ['C:\Users\Mason\anaconda3\_conda run -n ', conda_environment, ' ', 'python '  python_script, ' ', arguments];

    [status, result] = system(command);
    disp(result)

    

    folder = './diags/hdf5/';
    addpath(folder)

    files = dir(folder);
    files = files(~[files.isdir]);

    len = length(data);

    for ii = 1:numel(files)
        diag_beam = FBPIC_ReadElectron(folder, files(ii).name,'electrons_witness');
        lucBeam = YABP_FBPIC2YABP(diag_beam);
        lucBeam.BunchInterval = 1;
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







    






