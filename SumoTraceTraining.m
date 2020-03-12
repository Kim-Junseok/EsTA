%% Variables

numNodes = 155;
numSamples = 20;
validNum = 100000;
testNum = 100000;
scs = 60;

%% All data

for nodeInd = 1:numNodes
    fprintf(strcat('Get sample from node',num2str(nodeInd), " "));
    nodeData = rmmissing(readmatrix (strcat('data-', num2str(scs), '/node-', num2str(nodeInd), '.txt'))')';
    
    prevCellId = -1;
    cellSample = [];
    cellSampleSet = [];
    setInd = 1;
    
    for samInd = 1:length(nodeData)
        rssi = nodeData(samInd, 10);
        channelState = nodeData(samInd, 9);
        currCellId = nodeData(samInd, 11);
        taRegion = nodeData(samInd, 12);
        
        if prevCellId == currCellId
            temp = [rssi, channelState, taRegion];
            cellSample = [cellSample; temp];
        else
            prevCellId = currCellId;
            if ~isempty(cellSample)
                cellSampleSet(setInd).cellSample = cellSample;
                setInd = setInd + 1;
                cellSample = [];
            end
        end
    end
    
    cellSampleSet(setInd).cellSample = cellSample;
    sampleLength = 0;
    for ind = 1:length(cellSampleSet)
        sampleLength = sampleLength + (length(cellSampleSet(ind).cellSample)-numSamples+1);
    end
    fprintf(strcat(num2str(sampleLength), '\n'));
    
    tempDataSet = [];
    tempChannelLabelSet = [];
    tempTaLabelSet = [];

    fInd = 0;
    eInd = 0;
    for ind = 1:length(cellSampleSet)
        tempCellSample = cellSampleSet(ind).cellSample;    
        fInd = eInd;
        eInd = eInd + length(tempCellSample)-numSamples+1;
        for setInd = 1:length(tempCellSample)-numSamples+1
            tempDataSet(fInd+setInd,:) = tempCellSample(setInd:setInd+numSamples-1,1)';
            tempChannelLabelSet(fInd+setInd,:) = tempCellSample(setInd+numSamples-1,2);
            tempTaLabelSet(fInd+setInd,:) = tempCellSample(setInd+numSamples-1,3);
        end
    end
    
    if ~isempty(tempDataSet)
        if nodeInd == 1
            dataSet = tempDataSet;
            channelLabelSet = tempChannelLabelSet;
            taLabelSet = tempTaLabelSet;
        else
            dataSet = [dataSet; tempDataSet];
            channelLabelSet = [channelLabelSet; tempChannelLabelSet];
            taLabelSet = [taLabelSet; tempTaLabelSet];
        end
    end
end

fprintf(strcat('Sampling is completed, sample length: ', num2str(length(dataSet)), '\n'));

%% Divide data set in to train, valid, test sets

trainSet = dataSet(1:end-validNum-testNum,:);
trainChLabelSet = channelLabelSet(1:end-validNum-testNum,:);
trainTaLabelSet = taLabelSet(1:end-validNum-testNum,:);

validSet = dataSet(length(trainSet)+1:end-testNum,:);
validChLabelSet = channelLabelSet(length(trainSet)+1:end-testNum,:);
validTaLabelSet = taLabelSet(length(trainSet)+1:end-testNum,:);

testSet = dataSet(end-testNum+1:end,:);
testChLabelSet = channelLabelSet(end-testNum+1:end,:);
testTaLabelSet = taLabelSet(end-testNum+1:end,:);

%% LOS NLOS classification

trainSetLos = trainSet(trainChLabelSet==1,:);
trainSetNlos = trainSet(trainChLabelSet==0,:);
trainTaLabelSetLos = trainTaLabelSet(trainChLabelSet==1,:);
trainTaLabelSetNlos = trainTaLabelSet(trainChLabelSet==0,:);

fprintf('LOS \n')
for taInd=0:max(taLabelSet)+1
    list = find(trainTaLabelSetLos == taInd);
    if(~isempty(list))
        rsrpList = trainSetLos(list,end);
        fprintf(strcat(' TA', num2str(taInd), ','))
        figure(1)
        histfit(rsrpList(:));
        probCha.los(taInd+1) = fitdist(rsrpList(:), 'Normal');
        hold on
    end
end
fprintf('\n')

fprintf('NLOS \n')
for taInd=0:max(taLabelSet)+1
    list = find(trainTaLabelSetNlos == taInd);
    if(~isempty(list))
        rsrpList = trainSetNlos(list,end);
        fprintf(strcat(' TA', num2str(taInd), ','))
        figure(2)
        histfit(rsrpList(:));
        probCha.nlos(taInd+1) = fitdist(rsrpList(:), 'Normal');
        hold on
    end
end
fprintf('\n')


%% Print probability characteristics for all cells

fprintf('\nLOS RSRP mean values\n');
for ind = 1: length(probCha.los)
    fprintf(strcat('TA', num2str(ind),':',num2str(probCha.los(ind).mean),','));
end
fprintf('\n');
fprintf('NLOS RSRP mean values\n');
for ind = 1: length(probCha.nlos)
    fprintf(strcat('TA', num2str(ind),':',num2str(probCha.nlos(ind).mean),','));
end
fprintf('\n');

%%
% figure(1)
% for ind = 1:length(trainSet)
%     for taInd=0:max(taLabelSet)+1
%         list = find(trainTaLabelSetLos == taInd);
%         if(~isempty(list))
%             rsrpList = trainSetLos(list,end);
%             fprintf(strcat(' TA', num2str(taInd), ','))
%             figure(1)
%             histfit(rsrpList(:));
%             probCha.los(taInd+1) = fitdist(rsrpList(:), 'Normal');
%             hold on
%         end
%     end
%     plot(trainSet(ind,:), zeros(1,size(trainSet,2)), 'kx')
%     hold off
% end

%% Plotting depeding data dimension
% figure(3)
% for ind = 1:100
%     ind = randi(length(trainSetLos));
%     ta = taLabelSet(ind,:);
%     if ta == 0
%         plot(0, trainSetLos(ind, end), 'bo');
%     elseif ta == 1
%         plot(1, trainSetLos(ind, end), 'r+');
%     elseif ta == 2
%         plot(2, trainSetLos(ind, end), 'mx');
%     elseif ta == 3
%         plot(3, trainSetLos(ind, end), 'k*');
%     end
%     hold on
% end
% 
% figure(4)
% for ind = 1:100
%     ind = randi(length(trainSetLos));
%     ta = taLabelSet(ind,:);
%     if ta == 0
%         plot(trainSetLos(ind, end-1), trainSetLos(ind, end), 'bo');
%     elseif ta == 1
%         plot(trainSetLos(ind, end-1), trainSetLos(ind, end), 'r+');
%     elseif ta == 2
%         plot(trainSetLos(ind, end-1), trainSetLos(ind, end), 'mx');
%     elseif ta == 3
%         plot(trainSetLos(ind, end-1), trainSetLos(ind, end), 'k*');
%     end
%     hold on
% end
% 
% figure(5)
% for ind = 1:100
%     ind = randi(length(trainSetLos));
%     ta = taLabelSet(ind,:);
%     random = randi(size(trainSet, 2), 2, 1);
%     if ta == 0
%         plot3(trainSetLos(ind, random(1,:)), trainSetLos(ind, random(2,:)), trainSetLos(ind, end), 'bo');
%     elseif ta == 1
%         plot3(trainSetLos(ind, random(1,:)), trainSetLos(ind, random(2,:)), trainSetLos(ind, end), 'r+');
%     elseif ta == 2
%         plot3(trainSetLos(ind, random(1,:)), trainSetLos(ind, random(2,:)), trainSetLos(ind, end), 'mx');
%     elseif ta == 3
%         plot3(trainSetLos(ind, random(1,:)), trainSetLos(ind, random(2,:)), trainSetLos(ind, end), 'k*');
%     end
%     hold on
%     grid on
% end