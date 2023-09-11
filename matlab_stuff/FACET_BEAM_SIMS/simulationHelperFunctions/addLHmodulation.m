function [beamStructOut,ebeamout,deltagamma,t]=addLHmodulation(beamStructIn,Elaser,showplots)
% Hardcode FACET-II laser and undulator parameters
% Laser parameters
    lambda_laser = 760e-9;
    sigmar_laser = 200e-6;
    sigmat_laser = (9/2.355)*1e-12;
    Plaser = Elaser/sqrt(2*pi*sigmat_laser^2);
    offset = 0;% laser to e-beam offset if you want you can add it
% Undulator parameters
    K = 1.1699;
    lambdaw = 0.054;
    Nwig = 9;
% Electron beam
beamStructOut = beamStructIn;
ebeamIn = beamStructIn.Bunch.x;
%[x xp y yp z E (GeV) ] are coordinates of ebeamin
    x=ebeamIn(1,:)-mean(ebeamIn(1,:));    
    y=ebeamIn(3,:)-mean(ebeamIn(3,:));        
    gamma = ebeamIn(6,:)/0.511e-3;
    gamma0 = mean(ebeamIn(6,:))/0.511e-3;
    t =  (ebeamIn(5,:)-mean(ebeamIn(5,:)))/299792458;


lambda_r = lambdaw/2/gamma0^2*(1+K^2/2);% Assumes planar undulator
omega = 2*pi*299792458/lambda_r;% Resonant frequency
JJ = besselj(0,K^2/(4+2*K^2))-besselj(1,K^2/(4+2*K^2));
totalLaserEnergy = sqrt(2*pi*sigmat_laser^2)*Plaser;

% Laser is assumed Gaussian with peak power Plaser
% This formula from eq. 8 Huang PRSTAB 074401 2004
mod_amplitude = sqrt(Plaser/8.7e9)*K*lambdaw*Nwig/gamma0/sigmar_laser*JJ;
disp(num2str(mod_amplitude/sqrt(Plaser)))
%offset = 1.0;% temporal offset between laser and e-beam in units of laser wavelengths

% Calculate induced modulation deltagamma
deltagamma = mod_amplitude.*exp(-0.25.*(x.^2+y.^2)./sigmar_laser^2).*sin(omega.*t+offset*2*pi).*exp(-0.5.*((t-offset*sigmat_laser)./sigmat_laser).^2);    

ebeamout = ebeamIn;
ebeamout(6,:) = ebeamIn(6,:)/0.511e-3 + deltagamma; 

beamStructOut.Bunch.x(6,:) = 0.511e-3*ebeamout(6,:);

%% Plot stuff
if showplots
    z0 = 0;
    zmax = 1;
ind = ebeamIn(5,:)<(z0+zmax)*1e-3 & ebeamIn(5,:)>(z0-zmax)*1e-3;
figure(23)

subplot(3,2,1);
scatter(ebeamIn(5,:).*1e3,ebeamIn(6,:)*1e3,'.');box on;
hold on
plot(z0+zmax.*ones(100,1),linspace(min(ebeamIn(6,:)*1e3),max(ebeamIn(6,:)*1e3)),'--k')
plot(z0-zmax.*ones(100,1),linspace(min(ebeamIn(6,:)*1e3),max(ebeamIn(6,:)*1e3)),'--k')
xlabel('z [mm]');ylabel('Energy [MeV]');title('Before heater');

subplot(3,2,2);histogram(ebeamIn(6,ind).*1e3);box on;
xlabel('Energy [MeV]');ylabel('N');title(['std = ',num2str(std(ebeamIn(6,ind)*1e3),'%.1f'),' keV']);

subplot(3,2,3);scatter(ebeamout(5,:).*1e3,ebeamout(6,:)*0.511,'.');box on;
hold on
plot(z0+zmax.*ones(100,1),linspace(min(ebeamout(6,:)*0.511),max(ebeamout(6,:)*0.511)),'--k')
plot(z0-zmax.*ones(100,1),linspace(min(ebeamout(6,:)*0.511),max(ebeamout(6,:)*0.511)),'--k')
xlabel('z [mm]');ylabel('Energy [MeV]');title('After Heater');

subplot(3,2,4);histogram(ebeamout(6,ind).*0.511);box on;
xlabel('Energy [MeV]');ylabel('N');title(['std = ',num2str(std(ebeamout(6,ind)*0.511e3),'%.1f'),' keV']);

subplot(3,2,5);scatter(ebeamout(5,:).*1e3,deltagamma*0.511);box on;xlabel('z [mm]');ylabel('\delta E [MeV]');title('Modulation [zoom]');
xlim(z0+[-2,2].*lambda_r*1e3)

str = {['Und. K = ',num2str(K)],'',['Und. L [m] = ',num2str(Nwig*lambdaw)],'',['P_{laser} [MW]= ',num2str(Plaser.*1e-6)],...
    '',['E_{laser} [mJ]= ',num2str(sqrt(2*pi*sigmat_laser^2)*Plaser*1e3,'%.2f')],'',['\sigma_{r,laser} [um]= ',num2str(1e6*sigmar_laser,'%.0f')],...
    '',['FWHM_{t,laser} [ps]= ',num2str(1e12*2.355*sigmat_laser,'%.1f')]};
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','FontName'),'FontName','Times')
annotation('textbox',[0.6,0.1,0.14,0.2],'String',str,'FitBoxToText','on','FontSize',12)
%% Plot beam trace spaces
for i=1:6;psIn(:,i) = beamStructIn.Bunch.x(i,:);psOut(:,i) = beamStructOut.Bunch.x(i,:);end
%plot_beam_trace_spaces(psIn);% input beam
plot_beam_trace_spaces(psOut);% output beam
%% Plot the slice energy spread
% Remove particles that are stopped
nslices = 100;
dz = 0.1e-3;
for m=1:nslices
    zmin = min(ebeamout(5,:))+(m-1).*dz;
    zmax = zmin+dz;

    ind = ebeamout(5,:)<zmax & ebeamout(5,:)>zmin & ~beamStructIn.Bunch.stop;
    zval(m) = zmin;
    slice_espread_before(m) = std(ebeamIn(6,ind)*1e6);
    slice_espread_after(m) = std(ebeamout(6,ind)*0.511e3);
    npartRel(m) = sum(beamStructIn.Bunch.Q(ind));
end

figure;
yyaxis left
    plot(zval.*1e3,slice_espread_before,'LineWidth',2);grid on;hold on
    plot(zval.*1e3,slice_espread_after,'--','LineWidth',2);
    set(gca,'FontSize',20,'FontName','Times','LineWidth',2)   
    xlabel('z [mm]')
    ylabel('RMS Slice energy spread [keV]')
    
yyaxis right
    plot(zval.*1e3,npartRel*1e12,'LineWidth',2);    
    ylabel('Charge [pC]')
    xlim([-4,4])
    hold on
tvals = linspace(1.5*min(t),1.5*max(t),100);
laserProfile = exp(-0.5.*((tvals-offset*sigmat_laser)./sigmat_laser).^2);
    plot(tvals*3e8*1e3,100*laserProfile./max(laserProfile),'--')
    legend('Before Heater','After Heater','Charge','Laser Profile [a.u.]','Location','NorthWest')
      title(sprintf(['Laser E [mJ]= ',num2str(Elaser*1e3,'%.2f'),',  FWHM [ps] = ',num2str(1e12*sigmat_laser*2.355,'%.1f')]))
end