function [SimParams,SimStructs] = getWMMSEMatrix(SimParams,SimStructs)

switch SimParams.weightedSumRateMethod
    
    case 'PreScheduling'
        
        [SimParams,SimStructs] = getPreWeightedMMSEDesign(SimParams,SimStructs);

    case 'PerformScheduling'
        
        [SimParams,SimStructs] = getWeightedMMSEDesign(SimParams,SimStructs);
        
    case 'DistScheduling'
        
        [SimParams,SimStructs] = getDWeightedMMSEDesign(SimParams,SimStructs);
        
    case 'CNetworkBF'
        
        [SimParams,SimStructs] = getCNetworkBFWMMSEDesign(SimParams,SimStructs);
        
    case 'DNetworkBF'
        
        [SimParams,SimStructs] = getDNetworkBFWMMSEDesign(SimParams,SimStructs);
        
        
    case 'StreamScheduling'
        
        [SimParams,SimStructs] = getStrWeightedMMSEDesign(SimParams,SimStructs);
        
    otherwise
        
        display('Unknown Weighted Sum Rate Option !');
                
end

if strcmp(SimParams.DebugMode,'true')
    
    displayArray = zeros(1,SimParams.nBases + 1);
    
    for iBase = 1:SimParams.nBases
        displayArray(1,iBase) = trace(SimStructs.baseStruct{iBase,1}.P{1,1}' * SimStructs.baseStruct{iBase,1}.P{1,1});
    end
    
    displayArray(1,iBase + 1) = SimParams.sPower;
    display(displayArray);
    
end

end
