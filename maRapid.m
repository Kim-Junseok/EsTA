clear
close all

LOS
NLOS

% Set the Seed number
seed = 10;
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);
% End Setting Seed number

carrierFreq = 2; % 2 GHz
bandwidth = 20; % MHz
gnbTxPower = 30; % Tx power of gNB (dBm)
gnbAntGain = 17; % Antenna gain of gNB (dBi)
global gnbAntHeight; gnbAntHeight = 25; % gNB antenna height
ueNoiseFig = 9; % Noise figure (dB)
noisePsd = 10^((-174-30)*0.1) * 10^(ueNoiseFig*0.1);
gnbLocation = [0, 0, 25];
nRbs = 20; % number of resource blocks for SS block
cellRadius = 234;

iter = 1000;
%distUeGnb = zeros(1,iter);
%snrUeGnb = zeros(20,iter);

for seedInd = 1:seed
    for index=1:iter
        dist = 0;
        while dist < 35
            [x, y, z] = UeLocation(cellRadius*2); % UE location is (x, y, z)
            dist = CalcDist([x, y, z],  gnbLocation); % BS location is (0, 0, 25)
            uePosition.x = x;
            uePosition.y = y;
            uePosition.z = z;
        end
        
        twoDimDist = CalcDist([uePosition.x, uePosition.y], [0,0]);
        threeDimDist = dist;
        
        userInfo.position = uePosition;
        [los, bpDist] = DetermineLos(userInfo);
        
        if los == 1
            distUeGnbLos(:,iter*(seedInd-1)+index) = dist;
        else
            distUeGnbNlos(:,iter*(seedInd-1)+index) = dist;
        end
        
        %Tx power spectral density
        txPsd = 10^((gnbTxPower-30)*0.1) / (nRbs * 360000) * ones(nRbs, 1);
        
        % Antenna gain
        rxPsd = txPsd * 10^(gnbAntGain*0.1);
        
        % TR 36.873 V12.6.0 21ÆäÀÌÁö 3D-UMa Pathloss
        if los == 1  % LOS
            if twoDimDist <= bpDist
                pathLossDb = 22*log10(threeDimDist) + 28 + 20*log10(carrierFreq);
            else
                pathLossDb = 40*log10(threeDimDist) + 28 + 20*log10(carrierFreq) - 9*log10((bpDist)^2+(gnbLocation(3)-uePosition.z)^2);
            end
        else  % NLOS
            if twoDimDist <= bpDist
                pathLossDbLos = 22*log10(threeDimDist) + 28 + 20*log10(carrierFreq);
            else
                pathLossDbLos = 40*log10(threeDimDist) + 28 + 20*log10(carrierFreq) - 9*log10((bpDist)^2+(gnbLocation(3)-uePosition.z)^2);
            end
            pathLossDbNlos = 161.04 - 7.1*log10(20) + 7.5*log10(20) - (24.37 - 3.7*(20/gnbLocation(3))^2)*log10(gnbLocation(3)) + (43.42-3.1*log10(gnbLocation(3)))*(log10(threeDimDist)-3) + 20*log10(carrierFreq) - (3.2*(log10(17.625))^2-4.97) - 0.6*(uePosition.z - 1.5);
            
            pathLossDb = max(pathLossDbLos, pathLossDbNlos);
        end
        rxPsd = rxPsd .* 10^(-pathLossDb*0.1);
        
        % Shadowing
        if los == 1
            sfSigma = 4;
            sf = 10^(sfSigma*randn*0.1);
        else
            sfSigma = 6;
            sf = 10^(sfSigma*randn*0.1);
        end
        rxPsd = rxPsd * sf;
        
        if los == 1
            fastFading = losCh(16:35, index);
        else
            fastFading = nlosCh(16:35, index);
        end
        rxPsd = rxPsd .* 10.^(fastFading*0.1);
        
        noise = noisePsd * ones(20,1);
        snr  = rxPsd./noise;
        
        if los==1
            snrUeGnbLos(:,iter*(seedInd-1)+index) = 10*log10(snr);
        else
            snrUeGnbNlos(:,iter*(seedInd-1)+index) = 10*log10(snr);
        end
        
    end
end

distUeGnbLos = nonzeros(distUeGnbLos)';
distUeGnbNlos = nonzeros(distUeGnbNlos)';
snrUeGnbLos = nonzeros(mean(snrUeGnbLos))';
snrUeGnbNlos = nonzeros(mean(snrUeGnbNlos))';

figure(1)
plot(distUeGnbLos, snrUeGnbLos, 'b*')
hold on
plot(distUeGnbNlos(1:length(distUeGnbLos)), snrUeGnbNlos(1:length(distUeGnbLos)), 'r*')

len = length(distUeGnbLos);
taOneLos = [];
taTwoLos = [];
taThreeLos = [];
taOneNlos = [];
taTwoNlos = [];
taThreeNlos = [];

for index = 1:len
    if distUeGnbLos(index) <= 78
        taOneLos = [taOneLos, snrUeGnbLos(index)];
    elseif distUeGnbLos(index) <= 156
        taTwoLos = [taTwoLos, snrUeGnbLos(index)];
    else
        taThreeLos = [taThreeLos, snrUeGnbLos(index)];
    end
    
    if distUeGnbNlos(index) <= 78
        taOneNlos = [taOneNlos, snrUeGnbNlos(index)];
    elseif distUeGnbNlos(index) <= 156
        taTwoNlos = [taTwoNlos, snrUeGnbNlos(index)];
    else
        taThreeNlos = [taThreeNlos, snrUeGnbNlos(index)];
    end
end

taOne = [taOneLos taOneNlos];
taTwo = [taTwoLos taTwoNlos];
taThree = [taThreeLos taThreeNlos];

num = 10;

% figure(2)
% [p,x] = hist(taOneLos); plot(x,p/sum(p)); %PDF
% hold on
% [p,x] = hist(taTwoLos); plot(x,p/sum(p)); %PDF
% [p,x] = hist(taThreeLos); plot(x,p/sum(p)); %PDF

% figure(3)
% [p,x] = hist(taOne,num); plot(x,p/sum(p)); %PDF
% hold on
% [p,x] = hist(taTwo,num); plot(x,p/sum(p)); %PDF
% [p,x] = hist(taThree,num); plot(x,p/sum(p)); %PDF

figure(4)
histfit(taOneLos);
hold on
histfit(taTwoLos);
histfit(taThreeLos);

figure(5)
histfit(taOneNlos);
hold on
histfit(taTwoNlos);
histfit(taThreeNlos);

pdTaOneLos = fitdist(taOneLos','Normal');
pdTaTwoLos = fitdist(taTwoLos','Normal');
pdTaThreeLos = fitdist(taThreeLos','Normal');
pdTaOneNlos = fitdist(taOneNlos','Normal');
pdTaTwoNlos = fitdist(taTwoNlos','Normal');
pdTaThreeNlos = fitdist(taThreeNlos','Normal');

meanLos = [pdTaOneLos.mu pdTaTwoLos.mu pdTaThreeLos.mu]
stdLos = [pdTaOneLos.std pdTaTwoLos.std pdTaThreeLos.std]
meanNlos = [pdTaOneNlos.mu pdTaTwoNlos.mu pdTaThreeNlos.mu]
stdNlos = [pdTaOneNlos.std pdTaTwoNlos.std pdTaThreeNlos.std]

% figure(6)
% histfit(taOne, num, 'kernel');
% hold on
% histfit(taTwo, num, 'kernel');
% histfit(taThree, num, 'kernel');