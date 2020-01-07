period = 160; % SS block period (ms)
% velocity = 3; % velocity of UE (km/h)
cellRadius = 234;
velMaxInd = 2;
velocityList = [3 30];
samplesMaxInd = 10;
[X, Y] = meshgrid(velocityList, 5*(1:samplesMaxInd));

for velInd = 1:velMaxInd
    %velocity = 3*velInd;
    %veolocity = 3;
    velocity = velocityList(velInd);
    vel = velocity*1000/3600; % (m/s)
    
    iter = 10000;
    ueInfo = zeros(7, iter); % Location, LOS/NLOS, SNR of SS block
    
    % Set the Seed number
    seed = 3;
    s = RandStream('mt19937ar','Seed',seed);
    RandStream.setGlobalStream(s);
    % End Setting Seed number
    
    while dist < 35
        [x, y, z] = UeLocation(cellRadius*2); % UE location is (x, y, z)
        dist = CalcDist([x, y, z],  gnbLocation); % BS location is (0, 0, 25)
        uePosition.x = x;
        uePosition.y = y;
        uePosition.z = z;
    end
    
    % Make UE trace (location, snr of SS block, and timing advance value)
    ueInfo(3,:) = uePosition.z * ones(1, iter);
    trace = GaussRWconfinedCircleReflectingBoundaries(iter-1, 2000, vel * period*0.001, uePosition.x, uePosition.y);
    ueInfo(1:2,:) = trace';
    
    for index = 1:iter
        
        twoDimDist = CalcDist([ueInfo(1,index), ueInfo(2,index)], [0, 0]);
        threeDimDist = CalcDist([ueInfo(1,index), ueInfo(2,index), ueInfo(3,index)], gnbLocation);
        
        userInfo.position = uePosition;
        [los, bpDist] = DetermineLos(userInfo);
        
        %if los == 1
        %    distUeGnbLos(:,iter*(seedInd-1)+index) = dist;
        %else
        %    distUeGnbNlos(:,iter*(seedInd-1)+index) = dist;
        %end
        
        if los == 1
            ueInfo(4,index) = 1; % LOS
        else
            ueInfo(4,index) = 0; % NLOS
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
        
        ueInfo(5,index) = mean(10*log10(snr));
        
        % TA value
        if threeDimDist <= 78
            ueInfo(6,index) = 1;
        elseif threeDimDist <= 156
            ueInfo(6,index) = 2;
        else
            ueInfo(6,index) = 3;
        end
        
    end
    
    %% TA estimation
    
    % numSamples = 20; % The number of SNR samples for estimating timing advance value
    for numSamplesInd = 1:samplesMaxInd
        numSamples = numSamplesInd * 5;
        ueInfo(7,1:numSamples-1) = ueInfo(6,1:numSamples-1);
        
        for index = numSamples:iter
            temp = 0;
            for ind = 1:3 % timing advance value index
                likeliLos = normpdf(ueInfo(5,index-numSamples+1:index),meanLos(ind), stdLos(ind));
                if temp < prod(likeliLos)
                    temp = prod(likeliLos); % update temp probability
                    ueInfo(7,index) = ind; % update estimated timing advance value
                end
                
                likeliNlos = normpdf(ueInfo(5,index-numSamples+1:index),meanNlos(ind), stdNlos(ind));
                if temp < prod(likeliNlos)
                    temp = prod(likeliNlos); % update temp probability
                    ueInfo(7,index) = ind; % update estimated timing advance value
                end
            end
        end
        
        failSamples = length(nonzeros(ueInfo(7,:)-ueInfo(6,:)));
        totSamples = iter-numSamples+1;
        succProb(numSamplesInd, velInd) = (totSamples-failSamples)/totSamples*100;
        %(totSamples-failSamples)/totSamples*100;
    end
end

figure(2)
plot(Y(:,1), succProb(:,1), 'b')
hold on
plot(Y(:,2), succProb(:,2), 'r')
plot(Y(:,3), succProb(:,3), 'm')
% plot(Y(:,4), succProb(:,4))

figure(3)
surf(X,Y,succProb)
xlabel('velocity')
ylabel('Number of samples')