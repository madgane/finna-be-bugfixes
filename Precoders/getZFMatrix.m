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

        eP = pinv(augH);
        if SimParams.queueWt == 2
            [SimStructs.baseStruct{iBase}.P{iBand,1}] = performQueuedWF(eP,SimParams.sPower,Q);
        else
            [SimStructs.baseStruct{iBase}.P{iBand,1}] = performWFAlgorithm(eP,SimParams.sPower);
        end
        
    end
    
end
   
end
