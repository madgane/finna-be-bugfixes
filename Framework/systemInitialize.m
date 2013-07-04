

function [SimParams,SimStructs] = systemInitialize(SimParams,SimStructs)

for iUser = 1:SimParams.nUsers
    SimStructs.userStruct{iUser,1}.crThrpt = 1;
    SimStructs.userStruct{iUser,1}.PFmetric = 0;
    SimStructs.userStruct{iUser,1}.lastThrpt = 0;
    SimStructs.userStruct{iUser,1}.userID = iUser;
    SimStructs.userStruct{iUser,1}.tAllocation = 0;
    SimStructs.userStruct{iUser,1}.W = cell(SimParams.nBands,1);    
    SimStructs.userStruct{iUser,1}.dropThrpt = zeros(SimParams.nDrops,1);
    
    SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt = 0;
    SimStructs.userStruct{iUser,1}.trafficConfig.bufferLength = 'Inf';
    SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime = zeros(1,SimParams.nDrops);
end

for iBase = 1:SimParams.nBases
    SimStructs.baseStruct{iBase,1}.profile = 0;
    SimStructs.baseStruct{iBase,1}.baseID = iBase;
    SimStructs.baseStruct{iBase,1}.allocGains = cell(SimParams.nBands,1);
    SimStructs.baseStruct{iBase,1}.allocPattern = cell(SimParams.nBands,1);
    SimStructs.baseStruct{iBase,1}.P = cell(SimParams.nBands,1);
end

SimStructs.linkChan = cell(SimParams.nBases,SimParams.nBands);
SimStructs.actualChannel = cell(SimParams.nBases,SimParams.nBands);
SimStructs.chanHistory = cell(SimParams.nBases,SimParams.nBands);

if strcmp(SimParams.ChannelModel,'Jakes')
    for iUser = 1:SimParams.nUsers
        for iBase = 1:SimParams.nBases
            for iBand = 1:SimParams.nBands
                reset(SimStructs.JakesChStruct{iUser,iBase,iBand});
            end
        end
    end
end

SimParams.Debug.activeStatus = zeros(SimParams.nUsers,1);
SimParams.Debug.receivedRSSI = zeros(SimParams.nRxAntenna,SimParams.nRxAntenna,SimParams.nUsers);
end