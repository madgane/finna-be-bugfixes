
function [SimParams,SimStructs] = generateUserTrafficArrivals(SimParams,SimStructs)

enStatToolBox = 'true';

for iUser = 1:SimParams.nUsers
   
    SimStructs.userStruct{iUser,1}.trafficConfig.avgArrRate = SimParams.avgPktValues(1,iUser);
    
    cLambda = SimStructs.userStruct{iUser,1}.trafficConfig.avgArrRate;

    if strcmp(enStatToolBox,'true')
        poissonArrivals = [random('Poisson',cLambda,1,...
            (SimParams.nDrops - SimParams.gracePeriod)) zeros(1,SimParams.gracePeriod)];
    else
        poissonArrivals = [getPoisson(cLambda,1,...
            (SimParams.nDrops - SimParams.gracePeriod)) zeros(1,SimParams.gracePeriod)];
    end
    
    SimStructs.userStruct{iUser,1}.trafficHistory.pktArrival = poissonArrivals;
    
end

end