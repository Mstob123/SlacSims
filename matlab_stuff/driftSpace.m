function beamOut = driftSpace(beamIn, length)

            xOut = beamIn.x + beamIn.thtx*length;
            yOut = beamIn.y + beamIn.thty*length;
            zOut = beamIn.z + length;

            beamOut = beamInit(xOut, yOut, zOut, beamIn.thtx, beamIn.thty, beamIn.Pz, beamIn.gamma, beamIn.charge);
end