
clear all;

quadRange = linspace(-30,-36,10);

data = struct();

for ii=1:numel(quadRange)

    folder = ['fbOuts/fbOut' num2str(ii-1) '/hdf5/'];

    files = dir(folder);
    files = files(~[files.isdir]);

    data.mag(ii) = quadRange(ii);

    for jj = 1:numel(files)
        beam = FBPIC_ReadElectron(folder, files(jj).name,'electrons_witness');
        lucBeam = YABP_FBPIC2YABP(beam);
        data.beams(ii,jj) = lucBeam;
    end

end

[quads,diags] = size(data.beams);

%% Looped Plots

close all;


for ii=1:quads

    szxArray = zeros(1,diags);
    emitxArray = zeros(1,diags);
    szyArray = zeros(1,diags);
    emityArray = zeros(1,diags);
    dstArray = zeros(1,diags);
    

    for jj = 1:diags
        
        bm = data.beams(ii,jj);
        szxArray(jj) = std(bm.Bunch.x(1,:));
        szyArray(jj) = std(bm.Bunch.x(3,:));
        
        [nx,ny,nz]=GetNEmitFromBeam(bm,1);
        emitxArray(jj) = nx;
        emityArray(jj) = ny;

        dstArray(jj) = mean(bm.Bunch.x(5,:));
    end

    
    figure(1)
    p1(ii) = plot(dstArray,szxArray);
    hold on

    figure(2)
    p2(ii) = plot(dstArray,szyArray);
    hold on

    figure(3)
    p3(ii) = plot(dstArray,emitxArray);
    hold on

    figure(4)
    p4(ii) = plot(dstArray,emityArray);
    hold on

end

figure(1)
xlabel('z')
ylabel('\sigma_x')
legend(string(quadRange), "Location","northwest")

figure(2)
xlabel('z')
ylabel('\sigma_y')
legend(string(quadRange), "Location","northwest")

figure(3)
xlabel('z')
ylabel('\epsilon_x')
legend(string(quadRange), "Location","northwest")


figure(4)
xlabel('z [m]')
ylabel('\epsilon_y [rad m]')
legend(string(quadRange), "Location","northwest")





%% Other Plots



szxArray = zeros(1,quads);
emitxArray = zeros(1,quads);
szyArray = zeros(1,quads);
emityArray = zeros(1,quads);

emityiArray = zeros(1,quads);
emitxiArray = zeros(1,quads);

for ii = 1:quads
    bmi= data.beams(ii,1);
    szxArray(ii) = std(bmi.Bunch.x(1,:));
    szyArray(ii) = std(bmi.Bunch.x(3,:));

    bmf = data.beams(ii,end);
    [nxf,nyf,nzf]=GetNEmitFromBeam(bmf,1);
    emitxArray(ii) = nxf;
    emityArray(ii) = nyf;

    [nxi,nyi,nzi]=GetNEmitFromBeam(bmi,1);
    emitxiArray(ii) = nxi;
    emityiArray(ii) = nyi;


end

figure(5)
set(gcf, 'Position', [100, 100, 1600, 1000]);

subplot(2,2,2)
% plot(szxArray,emitxArray)
% hold on
% plot(szxArray,emityArray)
% plot(szxArray,sqrt(emitxArray.^2+emityArray.^2))
% xlabel('\sigma_x')


plot(quadRange,emitxArray)
hold on
plot(quadRange,emityArray)
plot(quadRange,sqrt(emitxArray.^2+emityArray.^2))

xlabel('quad strength')
ylabel('\epsilon')
legend('\epsilon_x','\epsilon_y','\epsilon_{total}')


subplot(2,2,3)
plot(szyArray,emitxArray)
hold on
plot(szyArray,emityArray)
plot(szyArray,sqrt(emitxArray.^2+emityArray.^2))
xlabel('\sigma_y')
ylabel('\epsilon')
legend('\epsilon_x','\epsilon_y','\epsilon_{total}')

subplot(2,2,4)
plot(sqrt(szxArray.^2+szyArray.^2),emitxArray)
hold on
plot(sqrt(szxArray.^2+szyArray.^2),emityArray)
plot(sqrt(szxArray.^2+szyArray.^2),sqrt(emitxArray.^2+emityArray.^2))

xlabel('\sigma_{total}')
ylabel('\epsilon')
legend('\epsilon_x','\epsilon_y','\epsilon_{total}')

subplot(2,2,1)
plot(quadRange,szxArray)
hold on
plot(quadRange,szyArray)
plot(quadRange,sqrt(szxArray.^2+szyArray.^2))

xlabel('quad strength')
ylabel('\sigma')
legend('\sigma_x','\sigma_y','\sigma_{total}')



figure(6)

beta0x = szxArray.^2./emitxiArray;
beta0y = szyArray.^2./emityiArray;

emitchngx = abs(emitxArray-emitxiArray);
emitchngy = abs(emityArray-emityiArray);


plot(quadRange, emitchngy, 'bo-')
xlabel('Quadrupole Strength [T/m]')
ylabel('\delta\epsilon_y [rad*m]')



%% tes

clear all;












    


    








