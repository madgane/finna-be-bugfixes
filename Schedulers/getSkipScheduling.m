function [SimParams,SimStructs] = getSkipScheduling(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
        
    for iBand = 1:SimParams.nBands        
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = uIndices;
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = ones(kUsers,1);
    end
end
end

