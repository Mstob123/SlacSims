function Luc2FBPICtxt(lucBeam, fileOut)

            bnch = lucBeam.Bunch.x;

            gam = bnch(6,:)*(1000/0.511);

            outMatrix = [bnch(1,:); bnch(3,:); -bnch(5,:); bnch(2,:).*gam; bnch(4,:).*gam; gam]';

            writematrix(outMatrix, fileOut, 'Delimiter', ' ')

end
