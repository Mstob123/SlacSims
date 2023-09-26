function ret = readData(data, plotSpots)

    for s = 1:length(data)
        
        z(s) = mean(data(s).beam.Bunch.x(5,:));
        [betax(s), betay(s)] = findBeta(data(s).beam);

        if (plotSpots)
            f = figure('Position',[100,100,600,600], 'Name', num2str(z(s)));
            beamImage(data(end).beam,[],[],[],[],[],[f, 0]);  
        end

    end

    f = figure('Position',[50,50,1200,800], 'Name', 'Beta');
    semilogy(z,betax,'b')
    hold on
    semilogy(z,betay,'r')
    legend('beta_x','beta_y')

    disp(z)


end