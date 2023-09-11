function [BEAMLINE]=changeLinacPhases(BEAMLINE,params)

% --- Set Linac Phases / degrees off-crest and amplitudes
idfbL2_1=findcells(BEAMLINE,'Name','K14_4*'); % don't touch E feedback stations
idfbL2_2=findcells(BEAMLINE,'Name','K14_5*'); 
idfbL3_1=[findcells(BEAMLINE,'Name','K19_1*') findcells(BEAMLINE,'Name','K19_2*') findcells(BEAMLINE,'Name','K19_3*')];
idfbL3_2=[findcells(BEAMLINE,'Name','K19_4*') findcells(BEAMLINE,'Name','K19_5*') findcells(BEAMLINE,'Name','K19_6*')];
for ip=[findcells(BEAMLINE,'Name','BEGL1F'):findcells(BEAMLINE,'Name','L1XFBEG') ...
    findcells(BEAMLINE,'Name','L1XFEND'):findcells(BEAMLINE,'Name','BEGBC11')]
  if isfield(BEAMLINE{ip},'Phase')         
    BEAMLINE{ip}.Phase=BEAMLINE{ip}.Phase+params.P1;   
    disp(['P1 = ',num2str(BEAMLINE{ip}.Phase)])
  end
    if isfield(BEAMLINE{ip},'Volt')      
    BEAMLINE{ip}.Volt=BEAMLINE{ip}.Volt*(1+params.V1);% Relative offset
  end
end
for ip=findcells(BEAMLINE,'Name','BEGL2F'):findcells(BEAMLINE,'Name','BEGBC14_1')
  if isfield(BEAMLINE{ip},'Phase') && ~ismember(ip,[idfbL2_1 idfbL2_2])% Change the linac phase for all but the E feedback stations
      
    BEAMLINE{ip}.Phase=BEAMLINE{ip}.Phase+params.P2;  
    disp(['P2 = ',num2str(BEAMLINE{ip}.Phase)])
  end
  if isfield(BEAMLINE{ip},'Volt') && ~ismember(ip,[idfbL2_1 idfbL2_2])% Change the linac phase for all but the E feedback stations
    BEAMLINE{ip}.Volt=BEAMLINE{ip}.Volt*(1+params.V2);% Relative offset  
    disp('changed L2 v for 2 bunch!')
  end
end
for ip=findcells(BEAMLINE,'Name','BEGL3F_1'):findcells(BEAMLINE,'Name','ENDL3F_1')
  if isfield(BEAMLINE{ip},'Phase') && ~ismember(ip,[idfbL3_1 idfbL3_2])% Change the linac phase for all but the E feedback stations
    BEAMLINE{ip}.Phase=0;% Don't scan L3 for now   
  end
  if isfield(BEAMLINE{ip},'Volt') && ~ismember(ip,[idfbL3_1 idfbL3_2])% Change the linac voltage for all but the E feedback stations
      %BEAMLINE{ip}.Volt=BEAMLINE{ip}.Volt*(1+4e-2);% Relative offset needed for 2 bunch case to get to 10 GeV
      disp(['V3 = ',num2str(BEAMLINE{ip}.Volt)])
      %disp('changed L3 v for 2 bunch!')
  end
end