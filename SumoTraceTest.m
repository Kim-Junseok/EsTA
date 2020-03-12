%% Input

threshold = 0;
esta = [];

%% Validation
for ind = 1:length(validSet)
    if validChLabelSet(ind,:) == 1 % LOS
        likelihood = inf(length(probCha.los),1); % Initialization for likelihood value for each TA region
        for taInd = 1:length(probCha.los)
            rsrpList = validSet(ind,1:end);
            probList = normpdf(rsrpList, probCha.los(taInd).mean, probCha.los(taInd).std);
            probList(probList < threshold) = []; % configurable
            if isempty(probList)
                likelihood(taInd) = -inf;
            else
                likelihood(taInd) = sum(log10(probList));
            end
        end
    else % NLOS
        likelihood = zeros(length(probCha.nlos),1); % Initialization for likelihood value for each TA region
        for taInd = 1:length(probCha.nlos)
            rsrpList = validSet(ind,1:end);
            probList = normpdf(rsrpList, probCha.nlos(taInd).mean, probCha.los(taInd).std);
            probList(probList < threshold) = []; % configurable
            if isempty(probList)
                likelihood(taInd) = -inf;
            else
                likelihood(taInd) = sum(log10(probList));
            end
        end
    end
    esta(ind,:) = find(likelihood==max(likelihood)) - 1;
end

%accuracy = sum(esta == validTaLabelSet)/length(validTaLabelSet);
accuracy = length(find(abs(esta - validTaLabelSet) <= 2))/length(validTaLabelSet);


