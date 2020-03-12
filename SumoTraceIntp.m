clear
S = load('berlin-trace-100.mat');
berlinTrace = S.berlintrace100;
len = length(berlinTrace);
nodeNum = 155;
for ind=1:nodeNum
    userTrace(nodeNum).node = [];
    userTrace(nodeNum).location = [];
    userTrace(nodeNum).velocity = [];
    userTraceInter(nodeNum).node = [];
    userTraceInter(nodeNum).location = [];
    userTraceInter(nodeNum).velocity = [];    
end

for ind=1:len
    nodeInd = berlinTrace(ind,2) + 1; % node Ind > 0
    if(isempty(userTrace(nodeInd).node)) % First
        userTrace(nodeInd).node = nodeInd;
        traceInd = 1;
        userTrace(nodeInd).location(traceInd,:) = [berlinTrace(ind, 1) berlinTrace(ind, 3:5)];
        userTrace(nodeInd).velocity(traceInd,:) = [berlinTrace(ind, 1) berlinTrace(ind, 6:8)];
    else
        traceInd = size(userTrace(nodeInd).location, 1)+ 1; % row size + 1
        userTrace(nodeInd).location(traceInd,:) = [berlinTrace(ind, 1) berlinTrace(ind, 3:5)];
        userTrace(nodeInd).velocity(traceInd,:) = [berlinTrace(ind, 1) berlinTrace(ind, 6:8)];
    end
end

figure(1)
for nodeInd=1:nodeNum
    plot(userTrace(nodeInd).location(:,2), userTrace(nodeInd).location(:,3), 'o')
    xlabel('X-axis (m)', 'FontSize', 11, 'FontName', 'Arial');
    ylabel('Y-axis (m)', 'FontSize', 11, 'FontName', 'Arial');
    hold on
end

% figure(2)
% for nodeInd=temp+1:nodeNum
%     plot(userTrace(nodeInd).location(:,2), userTrace(nodeInd).location(:,3), 'o')
%     xlabel('X-axis (m)', 'FontSize', 11, 'FontName', 'Arial');
%     ylabel('Y-axis (m)', 'FontSize', 11, 'FontName', 'Arial');
%     hold on
% end

%% Interpolation

step = 5; % 5 ms

for userInd=1:nodeNum
    originLoc = userTrace(userInd).location;
    originVel = userTrace(userInd).velocity;
    tempLoc = userTrace(userInd).location(1,:);
    tempVel = userTrace(userInd).velocity(1,:);
    
    timeInd = 1; %Index initialization
    originTimeInd = 1; %Index initialization
    time = originLoc(1,1); %Initial time
    while (tempLoc(timeInd,1) < originLoc(end,1))
        time = time + step * 1e6; % unit (ms->ns)
        if (time < originLoc(originTimeInd+1, 1))
            timeInd = timeInd + 1;
            tempLoc(timeInd, 1) = time;
            tempLoc(timeInd, 2:4) = tempLoc(timeInd-1, 2:4) + step * 1e-3 * originVel(originTimeInd, 2:4); % unit (ms->s)
            tempVel(timeInd, 1) = time;
            tempVel(timeInd, 2:4) = originVel(originTimeInd, 2:4);
        else
            timeInd = timeInd + 1;
            tempLoc(timeInd, 1) = time;
            tempLoc(timeInd, 2:4) = tempLoc(timeInd-1, 2:4) + step * 1e-3 * originVel(originTimeInd, 2:4); % unit (ms->s)
            originTimeInd = originTimeInd + 1;
            tempVel(timeInd, 1) = time;
            tempVel(timeInd, 2:4) = originVel(originTimeInd, 2:4);
        end
    end
    userTraceInter(userInd).node = userInd;
    userTraceInter(userInd).location = tempLoc;
    userTraceInter(userInd).velocity = tempVel;
end

figure(2)
for nodeInd=1:nodeNum
    plot(userTraceInter(nodeInd).location(:,2), userTraceInter(nodeInd).location(:,3), 'o')
    hold on
end