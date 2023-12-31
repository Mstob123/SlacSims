function data=emitDataUnpack(otruse,emitData)
data.energy=emitData(1);
data.emitx=emitData(2);
data.demitx=emitData(3);
data.emitxn=emitData(4);
data.demitxn=emitData(5);
data.embmx=emitData(6);
data.dembmx=emitData(7);
data.bmagx=emitData(8);
data.dbmagx=emitData(9);
data.bcosx=emitData(10);
data.dbcosx=emitData(11);
data.bsinx=emitData(12);
data.dbsinx=emitData(13);
data.betax=emitData(14);
data.dbetax=emitData(15);
data.bx0=emitData(16);
data.alphx=emitData(17);
data.dalphx=emitData(18);
data.ax0=emitData(19);
data.chi2x=emitData(20);
data.emity=emitData(21);
data.demity=emitData(22);
data.emityn=emitData(23);
data.demityn=emitData(24);
data.embmy=emitData(25);
data.dembmy=emitData(26);
data.bmagy=emitData(27);
data.dbmagy=emitData(28);
data.bcosy=emitData(29);
data.dbcosy=emitData(30);
data.bsiny=emitData(31);
data.dbsiny=emitData(32);
data.betay=emitData(33);
data.dbetay=emitData(34);
data.by0=emitData(35);
data.alphy=emitData(36);
data.dalphy=emitData(37);
data.ay0=emitData(38);
data.chi2y=emitData(39);
notr=0;
for iotr=find(otruse)
  notr=notr+1;
  data.ido(notr)=emitData(39+notr);
end
lid=emitData(39+notr+1);
for idata=1:lid
  data.id(idata)=emitData(39+notr+1+idata);
end
lS=emitData(39+notr+1+lid+1);
for idata=1:lS
  data.S(idata)=emitData(39+notr+1+lid+1+idata);
end
ind=39+notr+1+lid+1+lS+1;
data.sigxf=emitData(ind:ind+lid-1);ind=ind+lid;
data.sigyf=emitData(ind:ind+lid-1);ind=ind+lid;
data.DX=emitData(ind:ind+notr-1);ind=ind+notr;
data.dDX=emitData(ind:ind+notr-1);ind=ind+notr;
data.DPX=emitData(ind:ind+notr-1);ind=ind+notr;
data.dDPX=emitData(ind:ind+notr-1);ind=ind+notr;
data.DY=emitData(ind:ind+notr-1);ind=ind+notr;
data.dDY=emitData(ind:ind+notr-1);ind=ind+notr;
data.DPY=emitData(ind:ind+notr-1);ind=ind+notr;
data.dDPY=emitData(ind:ind+notr-1);ind=ind+notr;
data.dp=emitData(ind);ind=ind+1;
data.sigx=emitData(ind:ind+notr-1);ind=ind+notr;
data.dsigx=emitData(ind:ind+notr-1);ind=ind+notr;
data.sigy=emitData(ind:ind+notr-1);ind=ind+notr;
data.dsigy=emitData(ind:ind+notr-1);ind=ind+notr;
data.xf=emitData(ind:ind+notr-1);ind=ind+notr;
data.yf=emitData(ind:ind+notr-1);ind=ind+notr;
for n=1:notr
  data.R{n}=reshape(emitData(ind:ind+15),4,[]); % 4x4 matrices
  ind=ind+16;
end
data.IP.ex0=emitData(ind);ind=ind+1;
data.IP.bx0=emitData(ind);ind=ind+1;
data.IP.ax0=emitData(ind);ind=ind+1;
data.IP.ey0=emitData(ind);ind=ind+1;
data.IP.by0=emitData(ind);ind=ind+1;
data.IP.ay0=emitData(ind);ind=ind+1;
data.IP.sigx=emitData(ind);ind=ind+1;
data.IP.dsigx=emitData(ind);ind=ind+1;
data.IP.sigpx=emitData(ind);ind=ind+1;
data.IP.dsigpx=emitData(ind);ind=ind+1;
data.IP.betax=emitData(ind);ind=ind+1;
data.IP.dbetax=emitData(ind);ind=ind+1;
data.IP.alphx=emitData(ind);ind=ind+1;
data.IP.dalphx=emitData(ind);ind=ind+1;
data.IP.sigy=emitData(ind);ind=ind+1;
data.IP.dsigy=emitData(ind);ind=ind+1;
data.IP.sigpy=emitData(ind);ind=ind+1;
data.IP.dsigpy=emitData(ind);ind=ind+1;
data.IP.betay=emitData(ind);ind=ind+1;
data.IP.dbetay=emitData(ind);ind=ind+1;
data.IP.alphy=emitData(ind);ind=ind+1;
data.IP.dalphy=emitData(ind);ind=ind+1;
if (ind>length(emitData)),return,end
data.waist.Lx=emitData(ind);ind=ind+1;
data.waist.dLx=emitData(ind);ind=ind+1;
data.waist.betax=emitData(ind);ind=ind+1;
data.waist.dbetax=emitData(ind);ind=ind+1;
data.waist.sigx=emitData(ind);ind=ind+1;
data.waist.dsigx=emitData(ind);ind=ind+1;
data.waist.Ly=emitData(ind);ind=ind+1;
data.waist.dLy=emitData(ind);ind=ind+1;
data.waist.betay=emitData(ind);ind=ind+1;
data.waist.dbetay=emitData(ind);ind=ind+1;
data.waist.sigy=emitData(ind);ind=ind+1;
data.waist.dsigy=emitData(ind);ind=ind+1;
end