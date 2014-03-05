
function [SimParams,SimStructs] = generateUserTrafficArrivals(SimParams,SimStructs)

enStatToolBox = 'true';

queueModel = char(SimParams.arrivalDist);
uscore_index_q = find(queueModel == '_');

if ~isempty(uscore_index_q)
    queueType = queueModel(1:uscore_index_q(1,1) - 1);
    maxPktArrival = str2double(queueModel(uscore_index_q(1,1) + 1:end));
else
    queueType = queueModel;
    maxPktArrival = SimParams.maxArrival(1,SimParams.iPkt);
end

switch queueType
    case 'Uniform'
        randArrival = rand(1,SimParams.nUsers) * maxPktArrival;
    case 'Constant'
        randArrival = ones(1,SimParams.nUsers) * maxPktArrival;
    case 'Fixed'
        randArrival = SimParams.FixedPacketArrivals;
    case 'ConstFixed'
        randArrival = ones(1,SimParams.nUsers) * maxPktArrival;
    case 'ConstFlow'
        randArrival = SimParams.FixedPacketArrivals;
end

SimParams.avgPktValues = randArrival;

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
    
    if strcmp(queueType,'ConstFixed')
    poissonArrivals = SimParams.avgPktValues(1,iUser) * ones(1,length(poissonArrivals));
    end
    
    if strcmp(queueType,'ConstFlow')
    poissonArrivals = SimParams.avgPktValues(1,iUser) * ones(1,length(poissonArrivals));
    end
    
    SimStructs.userStruct{iUser,1}.trafficHistory.pktArrival = poissonArrivals;
    
end


end