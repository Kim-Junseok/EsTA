function [los, bpDist] = DetermineLos(userInfo)

global gnbAntHeight;

twoDimDist = userInfo.twoDimDist;
if twoDimDist <= 18
    %losProb = 1;
    los = 1;
else
    if userInfo.position.z <= 13
        cprime = 0;
    elseif userInfo.position.z <= 23
        cprime = power((userInfo.position.z - 13) / 10, 1.5);
    else
        fprintf('Fatal Error !!!!!, DetermineLos.m')
    end
    
    losProb = (18/twoDimDist + exp(-twoDimDist/63) * (1-18/twoDimDist)) * (1 + cprime * (5/4) * power(twoDimDist/100, 3) * exp(-twoDimDist/150));
    
    if rand > losProb
        los = 0; % NLOS
    else
        los = 1; % LOS
    end
end


%if los ==1
    if userInfo.position.z - 1.5 < 12
        effHeight = 1;
    else
        bndNum = (userInfo.position.z - 1.5 - 9) / 3;
        if bndNum <= 0
            fprintf('Error!!!!!, DetermineLos.m')
        end
        %effHeight = 3 * unidrnd(bndNum) + 9;
        effHeight = 3 * bndNum * rand + 9;
    end
    effGnbHeight = gnbAntHeight - effHeight;
    effUeHeight = userInfo.position.z - effHeight;
    bpDist = 4 * effGnbHeight * effUeHeight * (2e9/3e8);
%else
%    bpDist = NaN;
%end

end