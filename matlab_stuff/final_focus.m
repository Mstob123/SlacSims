function beamOut = final_focus(beamLucfile, ffparamsfile, qStr, outTxt)

    if ~exist('outTxt')
       outTxt = ['/home/mstobbe/FBPIC_Sims/fbIns/fbIn' num2str(qStr) '.txt']
    end

    mE_SI = 9.1094e-31;
    c0_SI = 299792458;
    qE_SI = 1.6022e-19;

    load(ffparamsfile, 'ffparams');
    beamLuc = load(beamLucfile);

    dLengths = ffparams(1,:);
    quadLengths = ffparams(2,:);
    quadStrengths = ffparams(3,:);

    beam = struct();

    beam.x = beamLuc.Bunch.x(1,:);
    beam.thtx  = beamLuc.Bunch.x(2,:);

    beam.y = beamLuc.Bunch.x(3,:);
    beam.thty  = beamLuc.Bunch.x(4,:);

    beam.z = -beamLuc.Bunch.x(5,:);

    E = beamLuc.Bunch.x(6,:);
    beam.gamma = E*1000/0.511;

    beam.charge = beamLuc.Bunch.Q;
    beam.Pz = beam.gamma*c0_SI*mE_SI;

    for ii = 1:6

        qLength = quadLengths(ii);
        qStrength = quadStrengths(ii);
        beam = quadrupole(beam, qStrength, qLength);
        

        dlength = dLengths(ii);
        beam = driftSpace(beam, dlength);
        
    end

    beamOut = beam;

    LucBeamOut = beamLuc;
    LucBeamOut.Bunch.x(1,:) = beam.x;
    LucBeamOut.Bunch.x(2,:) = beam.thtx;
    LucBeamOut.Bunch.x(3,:) = beam.y;
    LucBeamOut.Bunch.x(4,:) = beam.thty;
    LucBeamOut.Bunch.x(5,:) = -beam.z;

    convert.bstore.BEGINNING.Bunch = LucBeamOut.Bunch;
    save("/home/mstobbe/matlab_scripts/beamOut.mat", "-struct", "convert")

   
    Luc2FBPICtxt(LucBeamOut, outTxt)

    %save("outBeam.mat", "-struct", "LucBeamOut")

end