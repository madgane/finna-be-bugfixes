function [SimParams,SimStructs] = getZFMatrix(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    for iBand = 1:SimParams.nBands
    
        augH = [];
        Q = zeros(SimParams.muxRank,1);
        pickUsers = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1};
        pickStreams = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1};
        
        for iUser = 1:length(pickUsers)
            cUser = pickUsers(iUser,1);cStream = pickStreams(iUser,1);
            [W,~,~] = svd(SimStructs.linkChan{iBase,iBand}(:,:,cUser));
            augH = [augH ; W(:,cStream)' * SimStructs.linkChan{iBase,iBand}(:,:,cUser)];
            SimStructs.userStruct{cUser,1}.W{iBand,1}(:,cStream) = W(:,cStream);
            Q(iUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
        end

        eP = ((augH' * augH) + SimParams.N * eye(SimParams.nTxAntenna)) \ augH';
        switch SimParams.queueWt
            case 1
                [SimStructs.baseStruct{iBase}.P{iBand,1}] = performWFAlgorithm(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
            case 2
                [SimStructs.baseStruct{iBase}.P{iBand,1}] = performQueuedWF(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand),Q);
            otherwise
                 [SimStructs.baseStruct{iBase}.P{iBand,1}] = performWFAlgorithm(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
        end
        
    end
    
end
   
end
