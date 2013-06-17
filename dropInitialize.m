
function [SimParams,SimStructs] = dropInitialize(SimParams,SimStructs)

% Channel Generation + Addition of Path Loss

for iBand = 1:SimParams.nBands
    for iBase = 1:SimParams.nBases
        for iUser = 1:SimParams.nUsers
            
            PL = 10^(SimParams.PL_Profile(iBase,iUser) / 20);
            estError = complex(randn(SimParams.nRxAntenna,SimParams.nTxAntenna),randn(SimParams.nRxAntenna,SimParams.nTxAntenna)) * sqrt(0.5 * SimParams.estError);
           
            switch (SimParams.ChannelModel)
                
                case 'AWGN'
                    
                    SimStructs.actualChannel{iBase,iBand}(:,:,iUser) = complex(ones(SimParams.nRxAntenna,SimParams.nTxAntenna), ...
                        zeros(SimParams.nRxAntenna,SimParams.nTxAntenna)) * PL;
                    SimStructs.linkChan{iBase,iBand}(:,:,iUser) = sqrt(1 - SimParams.estError) * SimStructs.actualChannel{iBase,iBand}(:,:,iUser) + estError;
                                        
                case 'IID'
                    
                    SimParams.elapsedFBDuration(iUser,1) = 0;
                    SimStructs.actualChannel{iBase,iBand}(:,:,iUser) = complex(randn(SimParams.nRxAntenna,SimParams.nTxAntenna), ...
                        randn(SimParams.nRxAntenna,SimParams.nTxAntenna)) * sqrt(0.5) * PL;
                    
                    SimStructs.linkChan{iBase,iBand}(:,:,iUser) = sqrt(1 - SimParams.estError) * SimStructs.actualChannel{iBase,iBand}(:,:,iUser) + estError;
                    
                case 'Jakes'
                    
                    updateElapsedFeedbackDuration = 0;
                    [~,PathGains] = step(SimStructs.JakesChStruct{iUser,iBase,iBand},ones(SimParams.SFSymbols,SimParams.nTxAntenna));
                    SimStructs.actualChannel{iBase,iBand}(:,:,iUser) = squeeze(PathGains(SimParams.SFSymbols,end,:,:)).' * PL;
                    
                    if ~mod((SimParams.iDrop - 1),SimParams.updateFeedback(iUser,1))
                        SimStructs.linkChan{iBase,iBand}(:,:,iUser) = sqrt(1 - SimParams.estError) * SimStructs.actualChannel{iBase,iBand}(:,:,iUser) + estError;
                        if iBase * iBand == 1
                            updateElapsedFeedbackDuration = 1;
                        end
                    end
                    
                    if updateElapsedFeedbackDuration == 1
                        SimParams.elapsedFBDuration(iUser,1) = 0;
                    else
                        SimParams.elapsedFBDuration(iUser,1) = SimParams.elapsedFBDuration(iUser,1) + 1;
                    end                    
                    
            end
        end
    end
end

alphaPF = SimParams.PF_dur^-1;
for iUser = 1:SimParams.nUsers
    SimStructs.userStruct{iUser,1}.W = cell(SimParams.nBands,1);
    SimStructs.userStruct{iUser,1}.trafficConfig.currentArrival = SimStructs.userStruct{iUser,1}.trafficHistory.pktArrival(1,SimParams.iDrop);
    SimStructs.userStruct{iUser,1}.PFmetric = (1 - alphaPF) * SimStructs.userStruct{iUser,1}.crThrpt + alphaPF * SimStructs.userStruct{iUser,1}.lastThrpt;
end

for iBase = 1:SimParams.nBases
    SimStructs.baseStruct{iBase}.assignedUsers = cell(SimParams.nBands,1);
    SimStructs.baseStruct{iBase}.assignedStreams = cell(SimParams.nBands,1);
    SimStructs.baseStruct{iBase}.P = cell(SimParams.nBands,1);
end

% Traffic Related Calculations !

for iUser = 1:SimParams.nUsers
    currentBacklogs = SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt + SimStructs.userStruct{iUser,1}.trafficConfig.currentArrival;
    currentResidual =  currentBacklogs - floor(SimStructs.userStruct{iUser,1}.lastThrpt);
    
    SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt = max(0,currentResidual);
    SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime(1,SimParams.iDrop) = max(0,currentResidual);
    SimStructs.userStruct{iUser,1}.trafficHistory.pktService(SimParams.iPkt,SimParams.iDrop) = floor(SimStructs.userStruct{iUser,1}.lastThrpt);
    SimStructs.userStruct{iUser,1}.lastThrpt = 0;
end

activeUsers = zeros(SimParams.nUsers,1);
weighingFactor = zeros(SimParams.nUsers,1);

for iUser = 1:SimParams.nUsers
    weighingFactor(iUser,1) = SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt;
    activeUsers(iUser,1) = sign(SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt);
end

for iUser = 1:SimParams.nUsers
    if strcmp(SimParams.weighingEqual,'true')
        SimStructs.userStruct{iUser,1}.weighingFactor = activeUsers(iUser,1);
    else
        SimStructs.userStruct{iUser,1}.weighingFactor = weighingFactor(iUser,1);
    end
end

for iUser = 1:SimParams.nUsers
    switch SimParams.mdpFactor
        case 0
        case 2
            chnInstant = SimParams.elapsedFBDuration(iUser,1) * SimParams.sampTime;
            mdpFactor = abs(besselj(0,(2 * pi * SimParams.userDoppler(iUser,1) * chnInstant)));
            SimStructs.userStruct{iUser,1}.weighingFactor = SimStructs.userStruct{iUser,1}.weighingFactor * (mdpFactor)^(SimParams.elapsedFBDuration(iUser,1));
        otherwise
            SimStructs.userStruct{iUser,1}.weighingFactor = SimStructs.userStruct{iUser,1}.weighingFactor * (SimParams.mdpFactor)^(SimParams.elapsedFBDuration(iUser,1));
    end
end

end
