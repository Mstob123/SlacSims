function [quadGradient, Bvals,BDES] = getQuadBDES(BEAMLINE,istart,iend,plotFlag)
% Note about converting from quad gradient to BDES from Nick
%The gradient in T/m is given by c*BDES/L/energy/10
%Energy is given in eV
%BDES is in kG.m
%BDES [KG.m] = (10*L*energy/c) * quadGradient [T/m]; or equivalently:
%BDES [kG.m]= 2*10*Bvals[T]*energy/c; where Bvals is the Lucretia quadStruc.B
% According to Glen's scripts this is the conversion
%BDES = 2*10*BEAMLINE{ij}.B % The factor of 2 comes from the arrayfun part
quadEle = findcells(BEAMLINE,'Class','QUAD', istart, iend)
quadEle = quadEle(1:2:end)   ;
    for n=1:length(quadEle)
    Bvals(n) = BEAMLINE{quadEle(n)}.B;
    Bvals(n+1) = BEAMLINE{quadEle(n)}.B;
    svals(n) = BEAMLINE{quadEle(n)}.S;
    svals(n+1) = BEAMLINE{quadEle(n)}.S+BEAMLINE{quadEle(n)}.L;
    end
    % Convert Lucretia B vals to quad gradient in T/m and to BDES in kG.m
    energy = 125e6;% in the injector
    c = 3e8;
    quadGradient = Bvals/BEAMLINE{quadEle(n)}.L; % in T/m assumes all quads are the same length
    %BDES = 2*(10*BEAMLINE{quadEle(n)}.L*energy/c) * quadGradient;% Note the extra factor of 2 because all quads are the same length and are split in half In Lucretia
    %BDES = 10*Bvals*energy/c;% same as line above
    BDES = 2*10*Bvals;% This gets close (within a few percent) of Glen's values see getQuadBDES_GlensConversion.m in '/Users/cemma/Documents/Work/FACET-II/Lucretia_sims/E300/E300MeetingSims_062023'
if plotFlag
    figure;bar(svals,BDES);box on;grid on;
    xlabel('S [m]');ylabel('BDES [kG.m]');
    %figure;bar(svals,Bvals);box on;grid on;
    %xlabel('S [m]');ylabel('Lucretia B''L [T]');
    xlim([0,max(svals)])
    set(gca,'FontSize',20,'FontName','Times','LineWidth',2)
   
end
end