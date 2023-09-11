function Tw = plotSpotSizeFromAtoB(BEAMLINE,Initial,pointB,pointA)

Tw = TwissPlot(1,pointB,Initial,[-1 -1 0]);
In=TwissToInitial(Tw,pointA,Initial);
Tw = TwissPlot(pointA,pointB,In,[-1 -1 -1]);
    emitx = 4e-6;emity = 4e-6;
    CES = 3e-5;%Correlated energy spread at s10 injector
    sigmax = sqrt(Tw.betax.*emitx./(Tw.P./0.511e-3)+(Tw.etax*CES).^2);
    sigmay = sqrt(Tw.betay.*emity./(Tw.P./0.511e-3)+(Tw.etay*CES).^2);
    figure; plot(Tw.S,1e4*sigmax,'LineWidth',2);hold on;plot(Tw.S,1e4*sigmay,'LineWidth',2);ylabel('Value');
    plot(Tw.S,-1e2*Tw.etax,'LineWidth',2);hold on;plot(Tw.S,1e2*Tw.etay,'LineWidth',2);legend('\sigma_x [0.1 mm]','\sigma_y [0.1 mm]','-\eta_x [cm]','\eta_y [cm]','Location','NorthWest');xlabel('s [m]');set(gca,'FontSize',18,'LineWidth',2);
    xlim([Tw.S(1),Tw.S(end)]);grid on;
    
    sigmaxDispersion = (Tw.etax*CES).^2;
    sigmaBetax = Tw.betax.*emitx./(Tw.P./0.511e-3);
    figure; semilogy(Tw.S,sigmaxDispersion./sigmaBetax,'LineWidth',2);ylabel('\eta_x^2\sigma_\delta^2/\beta_x\epsilon_x');
    xlabel('s [m]');set(gca,'FontSize',18,'LineWidth',2);
    xlim([Tw.S(1),Tw.S(end)]);grid on;  
    %ylim([1,1e4])
end