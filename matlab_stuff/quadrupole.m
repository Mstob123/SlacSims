function beamOut = quadrupole(beamIn, dBydx, length)

            qE_SI = 1.6022e-19;
            
            zOut = beamIn.z + length;

            k = qE_SI*dBydx./beamIn.Pz;

            xOut = beamIn.x.*cos(sqrt(k).*length) + beamIn.thtx./sqrt(k).*sin(sqrt(k).*length);
            thtxOut = beamIn.x.*-sqrt(k).*sin(sqrt(k).*length)+beamIn.thtx.*cos(sqrt(k).*length);

            yOut = beamIn.y.*cosh(sqrt(k).*length) + beamIn.thty./sqrt(k).*sinh(sqrt(k).*length);
            thtyOut = beamIn.y.*sqrt(k).*sinh(sqrt(k).*length)+beamIn.thty.*cosh(sqrt(k).*length);

            beamOut = beamInit(xOut, yOut, zOut, thtxOut, thtyOut, beamIn.Pz, beamIn.gamma, beamIn.charge);
end