% function [SimParams,SimStructs] = getZFPrjMatrix(SimParams,SimStructs)

load base.mat;

for iBase = 1:SimParams.nBases
    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
    
    for iBand = 1:SimParams.nBands
        
        D = cell(kUsers,1);
        W = cell(kUsers,1);
        V = cell(kUsers,1);
        eG = zeros(SimParams.maxRank,kUsers);
        eH = SimStructs.linkChan(:,:,uIndices,iBase,iBand);
        for iUser = 1:kUsers
            d = svd(eH(:,:,iUser));
            eG(1:length(d),iUser) = d;
            [W{iUser,1}, D{iUser,1}, V{iUser,1}] = svd(eH(:,:,iUser));
            if (SimParams.nTxAntenna == 1)
                D{iUser,1} = D{iUser,1}(1,1);
            end
            if (SimParams.nRxAntenna == 1)
                D{iUser,1} = diag(D{iUser,1});
            end
        end
        
        Vm = [];
        for iUser = 1:kUsers
            Vm = [Vm V{iUser,1} * D{iUser,1}];
        end
                
        augH = zeros(SimParams.maxRank,SimParams.nTxAntenna);
        [~, sortIndex] = sort(eG(:),'descend');sortIndex = sortIndex - 1;
        
        for iRank = 1:SimParams.maxRank
            if (iRank > SimParams.muxRank)
                break;
            end
            cUser = floor(sortIndex(iRank,1) / SimParams.maxRank) + 1;
            cStream = mod(sortIndex(iRank,1),SimParams.maxRank) + 1;
            augH(iRank,:) = W{cUser,1}(:,cStream)' * eH(:,:,cUser);
            SimStructs.userStruct{uIndices(cUser,1)}.W(:,cStream,iBand) = W{cUser,1}(:,cStream);
        end
        
        eP = pinv(augH);
        SimStructs.baseStruct{iBase}.P(:,:,iBand) = performWFAlgorithm(eP,SimParams.sPower);
        SimStructs.baseStruct{iBase}.allocPattern(:,iBand) = sortIndex;
        SimStructs.baseStruct{iBase}.allocGains(:,:,iBand) = eG;
        
    end
    
end

