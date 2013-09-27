function [SimParams,SimStructs] = getBDScheduling(SimParams,SimStructs)

for iBase = 1:SimParams.nBases
    
    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
    
    schedUsers = zeros(min(SimParams.muxRank,kUsers),1);
    schedStreams = zeros(min(SimParams.muxRank,kUsers),1);
    
    charScheduling = char(SimParams.SchedType);
    uscore_index = find(charScheduling == '_');
    caseStudy = charScheduling(uscore_index + 1:end);
    
    for iBand = 1:SimParams.nBands
        
        eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);
        
        switch (caseStudy)
            
            case 'RNS'
                
                iIndex = 0;
                xLocs = zeros(kUsers * SimParams.maxRank,2);
                augE = [];
                
                for iUser = 1:kUsers
                    cUser = uIndices(iUser,1);
                    [U,~,~] = svd(eH(:,:,iUser));
                    if SimParams.queueWt
                        M = U' * eH(:,:,iUser) * (SimStructs.userStruct{cUser,1}.weighingFactor);
                    else
                        M = U' * eH(:,:,iUser) * sign(SimStructs.userStruct{cUser,1}.weighingFactor);
                    end
                    for iRank = 1:SimParams.maxRank
                        iIndex = iIndex + 1;
                        augE = [augE M(iRank,:).'];
                        xLocs(iIndex,:) = [cUser iRank];
                    end
                end
                
                G = [];X = augE;
                for iStream = 1:min(SimParams.muxRank,kUsers)
                    ppVolume = zeros(kUsers,1);
                    for iUser = 1:kUsers
                        if iStream == 1
                            ppVolume(iUser,1) = norm(X(:,iUser));
                        else
                            V = repmat(X(:,iUser),1,iStream - 1);
                            V = sqrt(diag(V' * V));
                            K = abs(G' * X(:,iUser));
                            ppVolume(iUser,1) = prod(V - K);
                        end
                    end
                    
                    [~,sortI] = sort(ppVolume,'descend');
                    G = [G X(:,sortI(1,1)) / norm(X(:,sortI(1,1)))];
                    
                    schedUsers(iStream,1) = xLocs(sortI(1,1),1);
                    schedStreams(iStream,1) = xLocs(sortI(1,1),2);
                    
                end
                
                SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = schedStreams;
                
            case 'SP'
                
                iIndex = 0;
                xLocs = zeros(kUsers * SimParams.maxRank,2);
                augE = [];
                
                for iUser = 1:kUsers
                    cUser = uIndices(iUser,1);
                    [U,~,~] = svd(eH(:,:,iUser));
                    %[U,~,~] = iterative_svd(eH(:,:,iUser),4);
                    if SimParams.queueWt
                        M = U' * eH(:,:,iUser) * (SimStructs.userStruct{cUser,1}.weighingFactor);
                    else
                        M = U' * eH(:,:,iUser) * sign(SimStructs.userStruct{cUser,1}.weighingFactor);
                    end
                    for iRank = 1:SimParams.maxRank
                        iIndex = iIndex + 1;
                        augE = [augE M(iRank,:).'];
                        xLocs(iIndex,:) = [cUser iRank];
                    end
                end
                
                [~,~,sortA] = qr(augE,0);
                for iRank = 1:min(SimParams.muxRank,kUsers)
                    schedUsers(iRank,1) = xLocs(sortA(1,iRank),1);
                    schedStreams(iRank,1) = xLocs(sortA(1,iRank),2);
                end
                
                SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = schedStreams;
                
                
            case 'SPSS'
                
                nStreams = min(SimParams.maxRank,SimParams.nRxAntenna);
                augE = [];
                
                for iUser = 1:kUsers
                    [U,~,~] = svd(eH(:,:,iUser));
                    if SimParams.queueWt
                        M = U' * eH(:,:,iUser);M  = M.' * (SimStructs.userStruct{uIndices(iUser,1),1}.weighingFactor);
                    else
                        M = U' * eH(:,:,iUser);M  = M.' * sign(SimStructs.userStruct{uIndices(iUser,1),1}.weighingFactor);
                    end
                    augE = [augE M(:,1:nStreams)];
                end
                
                [~,~,sortA] = qr(augE,'vector');
                sIndicesSort = mod((sortA - 1),nStreams) + 1;
                uIndicesSort = floor((sortA - 1)/nStreams) + 1;
                totalCnt = min(SimParams.muxRank,length(sortA));
                
                if strcmp(SimParams.weightedSumRateMethod,'StreamScheduling')
                    totalCnt = max(totalCnt,2 * SimParams.maxRank - 1);
                    totalCnt = min(totalCnt,length(sortA));
                end
                
                SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = uIndices(uIndicesSort(1,1:totalCnt),1);
                SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = sIndicesSort(1,1:totalCnt)';
                
            case 'M-SP'
                
                iIndex = 0;
                xLocs = zeros(kUsers * SimParams.maxRank,2);
                augE = [];
                
                for iUser = 1:kUsers
                    cUser = uIndices(iUser,1);
                    [U,~,~] = svd(eH(:,:,iUser));
                    if SimParams.queueWt
                        M = U' * eH(:,:,iUser) * (SimStructs.userStruct{cUser,1}.weighingFactor);
                    else
                        M = U' * eH(:,:,iUser) * sign(SimStructs.userStruct{cUser,1}.weighingFactor);
                    end
                    for iRank = 1:SimParams.maxRank
                        iIndex = iIndex + 1;
                        augE = [augE M(iRank,:).'];
                        xLocs(iIndex,:) = [cUser iRank];
                    end
                end
                
                G = [];X = augE;
                for iStream = 1:min(SimParams.muxRank,kUsers)
                    ppVolume = zeros(kUsers,1);
                    for iUser = 1:kUsers
                        if iStream == 1
                            ppVolume(iUser,1) = norm(X(:,iUser));
                        else
                            V = G * inv(G' * G) * G';
                            K = X(:,iUser)' * V' * X(:,iUser) / (norm(X(:,iUser)) * norm(V' * X(:,iUser)));
                            ppVolume(iUser,1) = abs(acos(K)) * norm((eye(SimParams.nTxAntenna) - G * inv(G' * G) * G')' * X(:,iUser));
                        end
                    end
                    
                    [~,sortI] = sort(ppVolume,'descend');
                    schedUsers(iStream,1) = xLocs(sortI(1,1),1);
                    schedStreams(iStream,1) = xLocs(sortI(1,1),2);
                    G = [G X(:,sortI(1,1))];
                end
                
                SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = schedStreams;
                
            case 'SPIF'
                
                iIndex = 0;
                xLocs = zeros(kUsers * SimParams.maxRank,2);
                augE = [];
                
                userGains = zeros(kUsers,1);
                for iUser = 1:kUsers
                    userGains(iUser,1) = log(trace(SimParams.Debug.receivedRSSI(:,:,uIndices(iUser,1),iBand)));
                    if userGains(iUser,1) == Inf
                        userGains(iUser,1) = 1;
                    end
                end
                
                nUsersToSelect = min(SimParams.iDrop,SimParams.muxRank);
                
                for iUser = 1:kUsers
                    
                    cUser = uIndices(iUser,1);
                    [U,~,~] = svd(eH(:,:,iUser));
                    
                    if SimParams.queueWt
                        M = U' * eH(:,:,iUser) * (SimStructs.userStruct{cUser,1}.weighingFactor);
                    else
                        M = U' * eH(:,:,iUser) * sign(SimStructs.userStruct{cUser,1}.weighingFactor);
                    end
                    
                    M = M * userGains(iUser,1);
                    
                    for iRank = 1:SimParams.maxRank
                        iIndex = iIndex + 1;
                        augE = [augE M(iRank,:).'];
                        xLocs(iIndex,:) = [cUser iRank];
                    end
                end
                
                [~,~,sortA] = qr(augE,0);
                for iRank = 1:min(nUsersToSelect,kUsers)
                    schedUsers(iRank,1) = xLocs(sortA(1,iRank),1);
                    schedStreams(iRank,1) = xLocs(sortA(1,iRank),2);
                end
                
                SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = schedStreams;
                
            case 'SPFP'
                
                iIndex = 0;
                xLocs = zeros(kUsers * SimParams.maxRank,2);
                augE = [];
                
                for iUser = 1:kUsers
                    cUser = uIndices(iUser,1);
                    [U,~,~] = iterative_svd_fp(eH(:,:,iUser),3,16,4);
                    if SimParams.queueWt
                        M = U' * eH(:,:,iUser) * (SimStructs.userStruct{cUser,1}.weighingFactor);
                    else
                        M = U' * eH(:,:,iUser) * sign(SimStructs.userStruct{cUser,1}.weighingFactor);
                    end
                    for iRank = 1:SimParams.maxRank
                        iIndex = iIndex + 1;
                        augE = [augE M(iRank,:).'];
                        xLocs(iIndex,:) = [cUser iRank];
                    end
                end
                
                [~,~,sortA] = qr(augE,0);
                for iRank = 1:min(SimParams.muxRank,kUsers)
                    schedUsers(iRank,1) = xLocs(sortA(1,iRank),1);
                    schedStreams(iRank,1) = xLocs(sortA(1,iRank),2);
                end
                
                SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = schedStreams;
                
                
        end
        
    end
    
end

end
