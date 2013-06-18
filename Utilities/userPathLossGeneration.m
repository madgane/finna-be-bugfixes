
function [SimParams, SimStructs] = userPathLossGeneration(SimParams,SimStructs)

losDetails = cell(SimParams.nUsers,length(SimParams.wrapExtCells));
rssiDetails = zeros(length(SimParams.wrapExtCells),SimParams.nUsers);

for iBase = 1:length(SimParams.wrapExtCells)
    for iUser = 1:SimParams.nUsers
        
        separationM = abs(SimParams.wrapExtCells(iBase,1) - SimStructs.userStruct{iUser,1}.loc);
        
        antGain = getAntennaPatterGain(SimParams.wrapExtCells(iBase,1),...
            SimStructs.userStruct{iUser,1}.loc,SimParams.sysConfig.layoutFeatures);
        [avgRxGain, isLOS] = evaluateLTE_PL(SimParams,separationM);
        
        if isLOS
            losDetails{iUser,iBase} = 'true';
        else
            losDetails{iUser,iBase} = 'false';
        end
        
        rssiDetails(iBase,iUser) = avgRxGain + antGain;
        
    end
end

for iUser = 1:SimParams.nUsers    
    [sortVal,sortInd] = sort(rssiDetails(:,iUser),'descend');
    for iBase = 1:SimParams.nBases
        cBase = mod(sortInd - 1,SimParams.nBases) + 1;
        SimParams.PL_Profile(cBase,iUser) = sortVal(iBase,1);
        SimStructs.userStruct{iUser,1}.losFading{cBase,1} = losDetails(iUser,sortInd);
    end
end
