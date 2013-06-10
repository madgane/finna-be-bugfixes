function [SimParams,SimStructs] = getNoScheduling(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
        
    for iBand = 1:SimParams.nBands
        
        assUsers = conv(upsample(uIndices,SimParams.maxRank),ones(SimParams.maxRank,1));
        SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = assUsers(1:(end - SimParams.maxRank + 1));
        SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = repmat((1:SimParams.maxRank)',kUsers,1);

    end
end
end
