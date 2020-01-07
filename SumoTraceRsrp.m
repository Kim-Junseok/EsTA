LOS
NLOS

% Set the Seed number
seed = 1;
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);
% End Setting Seed number

% Parameters
carrierFreq = 2; % 2 GHz
bandwidth = 20; % MHz
gnbTxPower = 30; % Tx power of gNB (dBm)
gnbAntGain = 17; % Antenna gain of gNB (dBi)
global gnbAntHeight; gnbAntHeight = 25; % gNB antenna height
ueNoiseFig = 9; % Noise figure (dB)
nRbs = 20; % number of resource blocks for SS block
noisePsd = 10^((-174-30)*0.1) * 10^(ueNoiseFig*0.1);
scs = 15; % kHz

switch(scs)
    case 15
        mu = 0;
    case 30
        mu = 1;
    case 60
        mu = 2;
    case 120
        mu = 3;
    case 240
        mu = 4;
    case 480
        mu = 5;
end

% Topology
isd = 500;
center = [750, 400];
gnbLocation = [[center, 25];
    center(1), center(2)+isd, 25;
    center(1), center(2)-isd, 25;
    center(1)-isd/2*sqrt(3), center(2)+isd/2, 25;
    center(1)-isd/2*sqrt(3), center(2)-isd/2, 25;
    center(1)+isd/2*sqrt(3), center(2)+isd/2, 25;
    center(1)+isd/2*sqrt(3), center(2)-isd/2, 25]; % Compare 7 cells



for userInd = 1:length(userTraceInter)
    chInd = randi(40000, length(gnbLocation), 1); % Channel index for each gNB
    for timeInd = 1:length(userTraceInter(userInd).location)
        chInd = chInd + 1;
        for ind = 1:length(chInd)
            if chInd(ind) > 40000
                chInd(ind) = mod(chInd(ind), 40000); % Maximum index of channel trace is 40000
            end
        end
        
        uePosition.x = userTraceInter(userInd).location(timeInd, 2);
        uePosition.y = userTraceInter(userInd).location(timeInd, 3);
        %uePosition.z = userTraceInter(userInd).location(timeInd, 4);
        uePosition.z = 1.5;
        userTraceInter(userInd).location(timeInd, 4) = uePosition.z;
        
        for gnbInd = 1:length(gnbLocation)
            threeDimDistPerGnbs(gnbInd) = CalcDist([uePosition.x, uePosition.y, uePosition.z],  gnbLocation(gnbInd,:));
            twoDimDistPerGnbs(gnbInd) = CalcDist([uePosition.x, uePosition.y], gnbLocation(gnbInd, 1:2));
            
            % User information (2D dist, 3D dist, UE position)
            userInfo.threeDimDist = threeDimDistPerGnbs(gnbInd);
            userInfo.twoDimDist = twoDimDistPerGnbs(gnbInd);
            userInfo.position = uePosition;
            [losPerGnbs(gnbInd), bpDist] = DetermineLos(userInfo);
            
            %Tx power spectral density
            txPsd = 10^((gnbTxPower-30)*0.1) / (nRbs * scs * 1000 * 12) * ones(nRbs, 1);
            
            % Antenna gain
            rxPsd = txPsd * 10^(gnbAntGain*0.1);
            
            % TR 36.873 V12.6.0 21∆‰¿Ã¡ˆ 3D-UMa Pathloss
            if losPerGnbs(gnbInd) == 1  % LOS
                if userInfo.twoDimDist <= bpDist
                    pathLossDb = 22*log10(userInfo.threeDimDist) + 28 + 20*log10(carrierFreq);
                else
                    pathLossDb = 40*log10(userInfo.threeDimDist) + 28 + 20*log10(carrierFreq) - 9*log10((bpDist)^2+(gnbLocation(gnbInd, 3)-uePosition.z)^2);
                end
            else  % NLOS
                if userInfo.twoDimDist <= bpDist
                    pathLossDbLos = 22*log10(userInfo.threeDimDist) + 28 + 20*log10(carrierFreq);
                else
                    pathLossDbLos = 40*log10(userInfo.threeDimDist) + 28 + 20*log10(carrierFreq) - 9*log10((bpDist)^2+(gnbLocation(gnbInd, 3)-uePosition.z)^2);
                end
                pathLossDbNlos = 161.04 - 7.1*log10(20) + 7.5*log10(20) - (24.37 - 3.7*(20/gnbLocation(gnbInd, 3))^2)*log10(gnbLocation(gnbInd, 3)) + (43.42-3.1*log10(gnbLocation(gnbInd, 3)))*(log10(userInfo.threeDimDist)-3) + 20*log10(carrierFreq) - (3.2*(log10(17.625))^2-4.97) - 0.6*(uePosition.z - 1.5);
                
                pathLossDb = max(pathLossDbLos, pathLossDbNlos);
            end
            rxPsd = rxPsd .* 10^(-pathLossDb*0.1);
            
            % Shadowing
            if losPerGnbs(gnbInd) == 1
                sfSigma = 4;
                sf = 10^(sfSigma*randn*0.1);
            else
                sfSigma = 6;
                sf = 10^(sfSigma*randn*0.1);
            end
            rxPsd = rxPsd * sf;
            
            if losPerGnbs(gnbInd) == 1
                fastFading = losCh(16:35, chInd(gnbInd));
            else
                fastFading = nlosCh(16:35, chInd(gnbInd));
            end
            rxPsd = rxPsd .* 10.^(fastFading*0.1);
            
            % convert PSD [W/Hz] to linear power [W] for the single RE
            powerTxW = rxPsd * scs;
            rsrpDbm = 10 * log10(1000 * sum(powerTxW) / nRbs);
            
            rsrpDbmPerGnbs(gnbInd) = rsrpDbm;
        end
        
        %connecGnbInd = find(rsrpDbmPerGnbs==max(rsrpDbmPerGnbs)); %Connect with highest signal strength gNB
        connecGnbInd = find(threeDimDistPerGnbs==min(threeDimDistPerGnbs)); % Coneect with cloest gNB
        userTraceInter(userInd).location(timeInd, 5:7) = userTraceInter(userInd).velocity(timeInd, 2:4);
        userTraceInter(userInd).location(timeInd, 8) = threeDimDistPerGnbs(connecGnbInd);
        userTraceInter(userInd).location(timeInd, 9) = losPerGnbs(connecGnbInd);
        userTraceInter(userInd).location(timeInd, 10) = rsrpDbmPerGnbs(connecGnbInd);
        userTraceInter(userInd).location(timeInd, 11) = connecGnbInd;    
        
        % Timing advance
        taGranu = 16 * 64 / 2^(mu);
        unitDist = taGranu / (480 * 10^3 * 4096) * 10^6 * 3 * 10^2 / 2;
        taRegion = floor(threeDimDistPerGnbs(connecGnbInd) / unitDist);
        
        userTraceInter(userInd).location(timeInd, 12) = taRegion;
        
    end
    filename = strcat ('node-', num2str(userInd),'.txt');
    fileID = fopen(filename,'w');
    nbytes = fprintf(fileID,'%f %f %f %f %f %f %f %f %d %f %d %d \n',userTraceInter(userInd).location');
    fclose(fileID);
end


% distUeGnbLos = nonzeros(distUeGnbLos)';
% distUeGnbNlos = nonzeros(distUeGnbNlos)';
% snrUeGnbLos = nonzeros(mean(snrUeGnbLos))';
% snrUeGnbNlos = nonzeros(mean(snrUeGnbNlos))';
%
% figure(1)
% plot(distUeGnbLos, snrUeGnbLos, 'b*')
% hold on
% plot(distUeGnbNlos(1:length(distUeGnbLos)), snrUeGnbNlos(1:length(distUeGnbLos)), 'r*')
%
% len = length(distUeGnbLos);
% taOneLos = [];
% taTwoLos = [];
% taThreeLos = [];
% taOneNlos = [];
% taTwoNlos = [];
% taThreeNlos = [];
%
% for index = 1:len
%     if distUeGnbLos(index) <= 78
%         taOneLos = [taOneLos, snrUeGnbLos(index)];
%     elseif distUeGnbLos(index) <= 156
%         taTwoLos = [taTwoLos, snrUeGnbLos(index)];
%     else
%         taThreeLos = [taThreeLos, snrUeGnbLos(index)];
%     end
%
%     if distUeGnbNlos(index) <= 78
%         taOneNlos = [taOneNlos, snrUeGnbNlos(index)];
%     elseif distUeGnbNlos(index) <= 156
%         taTwoNlos = [taTwoNlos, snrUeGnbNlos(index)];
%     else
%         taThreeNlos = [taThreeNlos, snrUeGnbNlos(index)];
%     end
% end
%
% taOne = [taOneLos taOneNlos];
% taTwo = [taTwoLos taTwoNlos];
% taThree = [taThreeLos taThreeNlos];
%
% num = 10;
%
% % figure(2)
% % [p,x] = hist(taOneLos); plot(x,p/sum(p)); %PDF
% % hold on
% % [p,x] = hist(taTwoLos); plot(x,p/sum(p)); %PDF
% % [p,x] = hist(taThreeLos); plot(x,p/sum(p)); %PDF
%
% % figure(3)
% % [p,x] = hist(taOne,num); plot(x,p/sum(p)); %PDF
% % hold on
% % [p,x] = hist(taTwo,num); plot(x,p/sum(p)); %PDF
% % [p,x] = hist(taThree,num); plot(x,p/sum(p)); %PDF
%
% figure(4)
% histfit(taOneLos);
% hold on
% histfit(taTwoLos);
% histfit(taThreeLos);
%
% figure(5)
% histfit(taOneNlos);
% hold on
% histfit(taTwoNlos);
% histfit(taThreeNlos);
%
% pdTaOneLos = fitdist(taOneLos','Normal');
% pdTaTwoLos = fitdist(taTwoLos','Normal');
% pdTaThreeLos = fitdist(taThreeLos','Normal');
% pdTaOneNlos = fitdist(taOneNlos','Normal');
% pdTaTwoNlos = fitdist(taTwoNlos','Normal');
% pdTaThreeNlos = fitdist(taThreeNlos','Normal');
%
% meanLos = [pdTaOneLos.mu pdTaTwoLos.mu pdTaThreeLos.mu]
% stdLos = [pdTaOneLos.std pdTaTwoLos.std pdTaThreeLos.std]
% meanNlos = [pdTaOneNlos.mu pdTaTwoNlos.mu pdTaThreeNlos.mu]
% stdNlos = [pdTaOneNlos.std pdTaTwoNlos.std pdTaThreeNlos.std]
%
% % figure(6)
% % histfit(taOne, num, 'kernel');
% % hold on
% % histfit(taTwo, num, 'kernel');
% % histfit(taThree, num, 'kernel');