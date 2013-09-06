
function [SimParams,SimStructs] = getQWtdWSRMDist(SimParams,SimStructs)

cH = SimStructs.linkChan;
nBases = SimParams.nBases;
nBands = SimParams.nBands;

updatePrecoders = 'true';
usersPerCell = zeros(nBases,1);
cellUserIndices = cell(nBases,1);

% Debug Buffers initialization

SimParams.Debug.tempResource{2,1} = cell(SimParams.nUsers,1);
SimParams.Debug.tempResource{3,1} = cell(SimParams.nUsers,1);
SimParams.Debug.tempResource{4,1} = cell(SimParams.nUsers,SimParams.nBands);

for iBase = 1:nBases
    for iBand = 1:nBands
        cellUserIndices{iBase,1} = [cellUserIndices{iBase,1} ; SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}];
    end
    cellUserIndices{iBase,1} = unique(cellUserIndices{iBase,1});
    usersPerCell(iBase,1) = length(cellUserIndices{iBase,1});
end

nUsers = sum(usersPerCell);
QueuedPkts = zeros(nUsers,1);

for iBase = 1:nBases
    for iUser = 1:usersPerCell(iBase,1)
        cUser = cellUserIndices{iBase,1}(iUser,1);
        QueuedPkts(cUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
    end
end

switch SimParams.weightedSumRateMethod
    
    case 'PrimalDecomposition'
        
        alpha = 0.5;
        reIterate = 1;
        cvx_hist = -500 * ones(2,1);
        p_o = cell(nBases,1);q_o = cell(nBases,1);b_o = cell(nBases,1);cM = cell(nBases,1);
        
        xIndex = 0;
        for iBase = 1:nBases
            p_o{iBase,1} = ones(usersPerCell(iBase,1),nBands) / usersPerCell(iBase,1);
            q_o{iBase,1} = ones(usersPerCell(iBase,1),nBands) / usersPerCell(iBase,1);
            b_o{iBase,1} = ones(usersPerCell(iBase,1),nBands) * 10 + rand(usersPerCell(iBase,1),nBands);
        end
        
        fixedIF = ones(nUsers,nBases,nBands);
        
        while reIterate
            
            masterIterate = 1;totalCVXvalH = -50;
            
            while masterIterate
                
                cD = cell(nBases,1);
                
                for iBase = 1:nBases
                    
                    kUsers = usersPerCell(iBase,1);
                    
                    cvx_begin
                    
                    dual variables D{nUsers,nBands}
                    
                    expressions p(kUsers,nBands) q(kUsers,nBands)
                    variable M(SimParams.nTxAntenna,kUsers,nBands) complex
                    variables t(kUsers,nBands) b(kUsers,nBands) g(kUsers,nBands)
                    variables userObjective(kUsers,1) epiObjective
                    
                    minimize(epiObjective)
                    
                    subject to
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        abs(QueuedPkts(cUser,1) - sum(t(iUser,:))) <= userObjective(iUser,1);
                    end
                    
                    epiObjective >= norm(userObjective,1);
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            
                            intVector = sqrt(SimParams.N);
                            H = cH{iBase,iBand}(:,:,cUser);
                            
                            for jBase = 1:nBases
                                if jBase ~= iBase
                                    intVector = [intVector ; fixedIF(cUser,jBase,iBand)];
                                end
                            end
                            
                            for jUser = 1:kUsers
                                if jUser ~= iUser
                                    intVector = [intVector ; H * M(:,jUser,iBand)];
                                end
                            end
                            
                            D{cUser,iBand} : norm(intVector,2) <= sqrt(b(iUser,iBand));
                            log(1 + g(iUser,iBand)) >= t(iUser,iBand);
                            
                            p(iUser,iBand) = real(H * M(:,iUser,iBand));
                            q(iUser,iBand) = imag(H * M(:,iUser,iBand));
                            
                            q(iUser,iBand) == 0;
                            
                            (p_o{iBase,1}(iUser,iBand)^2 + q_o{iBase,1}(iUser,iBand)^2) / (b_o{iBase,1}(iUser,iBand)) + ...
                                (2 / b_o{iBase,1}(iUser,iBand)) * (p_o{iBase,1}(iUser,iBand) * (p(iUser,iBand) - p_o{iBase,1}(iUser,iBand))) + ...
                                (2 / b_o{iBase,1}(iUser,iBand)) * (q_o{iBase,1}(iUser,iBand) * (q(iUser,iBand) - q_o{iBase,1}(iUser,iBand))) - ...
                                (p_o{iBase,1}(iUser,iBand)^2 + q_o{iBase,1}(iUser,iBand)^2) / (2 * b_o{iBase,1}(iUser,iBand)^2) * ...
                                (b(iUser,iBand) - b_o{iBase,1}(iUser,iBand)) >= g(iUser,iBand);
                            
                        end
                        
                        norm(vec(M(:,:,iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                        
                        for iUser = 1:nUsers
                            intVector = [];
                            if SimStructs.userStruct{iUser,1}.baseNode ~= iBase
                                for jUser = 1:kUsers
                                    intVector = [intVector ; cH{iBase,iBand}(:,:,iUser) * M(:,jUser,iBand)];
                                end
                                D{iUser,iBand} : norm(intVector,2) <= fixedIF(iUser,iBase,iBand);
                            end
                        end
                        
                    end
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        norm(t(iUser,:),1) <= QueuedPkts(cUser,1) * log(2);
                    end
                    
                    cvx_end
                    
                    if strfind(cvx_status,'Solved')
                        cM{iBase,1} = M;
                        b_o{iBase,1} = full(b);
                        
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(iUser,:))) * log2(exp(1)),0);
                            SimParams.Debug.tempResource{2,1}{cUser,1} = [SimParams.Debug.tempResource{2,1}{cUser,1} sum(vec(t(iUser,:)))];
                            SimParams.Debug.tempResource{3,1}{cUser,1} = [SimParams.Debug.tempResource{3,1}{cUser,1} qDeviation];
                            for iBand = 1:nBands
                                SimParams.Debug.tempResource{4,1}{cUser,iBand} = [SimParams.Debug.tempResource{4,1}{cUser,iBand} t(iUser,iBand)];
                            end
                        end
                        
                    else
                        b_o{iBase,1} = b_o{iBase,1} * 2;
                    end
                    
                    cD{iBase,1} = D;
                    
                end
                
                fixedIFH = fixedIF;
                fixedIF = zeros(nUsers,nBases,nBands);
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:nUsers
                            baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                            if baseNode ~= iBase
                                fixedIF(iUser,iBase,iBand) = fixedIFH(iUser,iBase,iBand) + alpha * (cD{baseNode,1}{iUser,iBand} - cD{iBase,1}{iUser,iBand});
                            end
                        end
                    end
                end
                
                display(fixedIF);
                totalCVXval = 0;
                for iBase = 1:nBases
                    for iUser = 1:nUsers
                        for iBand = 1:nBands
                            totalCVXval = totalCVXval + max(cD{iBase,1}{iUser,iBand},0);
                        end
                    end
                end
                
                fixedIF = max(fixedIF,0);
                if abs(totalCVXvalH - totalCVXval) <= 1e-3
                    masterIterate = 0;
                else
                    totalCVXvalH = totalCVXval;
                end
                
            end
            
            for iBand = 1:nBands
                for iBase = 1:nBases
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        currentH = cH{iBase,iBand}(:,:,cUser);
                        p_o{iBase,1}(iUser,iBand) = real(currentH * cM{iBase,1}(:,iUser,iBand));
                        q_o{iBase,1}(iUser,iBand) = imag(currentH * cM{iBase,1}(:,iUser,iBand));
                    end
                end
            end
            
            if min(abs(cvx_optval - cvx_hist)) <= 1e-4
                reIterate = 0;
            else
                xIndex = xIndex + 1;
                cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
            end
            
        end
        
        updatePrecoders = 'false';
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = cM{iBase,1}(:,:,iBand);
            end
        end
        
    case 'PriDGenAlloc'
        
        alpha = 0.5;
        vW = cell(nUsers,nBands);
        maxRank = SimParams.maxRank;
        b_o = cell(nBases,1);cM = cell(nBases,1);
        p_o = cell(nBases,1);q_o = cell(nBases,1);
        
        reIterate = 1;
        cvx_hist = -500 * ones(2,1);
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    vW{cUser,iBand} = ones(SimParams.nRxAntenna,maxRank) / sqrt(SimParams.nRxAntenna);
                end
            end
            p_o{iBase,1} = ones(maxRank,usersPerCell(iBase,1),nBands) / usersPerCell(iBase,1);
            q_o{iBase,1} = ones(maxRank,usersPerCell(iBase,1),nBands) / usersPerCell(iBase,1);
            b_o{iBase,1} = ones(maxRank,usersPerCell(iBase,1),nBands) * 10 + rand(maxRank,usersPerCell(iBase,1),nBands);
        end
        
        xIndex = 0;
        fixedIF = ones(maxRank,nUsers,nBases,nBands);
        
        while reIterate
            
            masterIterate = 1;totalCVXvalH = -50;
            
            while masterIterate
                
                cD = cell(nBases,1);
                
                for iBase = 1:nBases
                    
                    kUsers = usersPerCell(iBase,1);
                    
                    cvx_begin
                    
                    dual variables D{maxRank,nUsers,nBands}
                    
                    expressions p(maxRank,kUsers,nBands) q(maxRank,kUsers,nBands)
                    variable M(SimParams.nTxAntenna,maxRank,kUsers,nBands) complex
                    variables t(maxRank,kUsers,nBands) b(maxRank,kUsers,nBands) g(maxRank,kUsers,nBands)
                    variables userObjective(kUsers,1) epiObjective
                    
                    minimize(epiObjective)
                    
                    subject to
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        abs(QueuedPkts(cUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1);
                    end
                    
                    epiObjective >= norm(userObjective,1);
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            
                            for iLayer = 1:maxRank
                                
                                intVector = sqrt(SimParams.N);
                                H = cH{iBase,iBand}(:,:,cUser);
                                
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        intVector = [intVector ; fixedIF(iLayer,cUser,jBase,iBand)];
                                    end
                                end
                                
                                for jUser = 1:kUsers
                                    if jUser ~= iUser
                                        for jLayer = 1:maxRank
                                            intVector = [intVector ; vW{cUser,1}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                        end
                                    else
                                        for jLayer = 1:maxRank
                                            if iLayer ~= jLayer
                                                intVector = [intVector ; vW{cUser,1}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                            end
                                        end
                                    end
                                end
                                
                                D{iLayer,cUser,iBand} : norm(intVector,2) <= sqrt(b(iLayer,iUser,iBand));
                                log(1 + g(iLayer,iUser,iBand)) >= t(iLayer,iUser,iBand);
                                
                                p(iLayer,iUser,iBand) = real(vW{cUser,1}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                q(iLayer,iUser,iBand) = imag(vW{cUser,1}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                
                                q(iLayer,iUser,iBand) == 0;
                                
                                (p_o{iBase,1}(iLayer,iUser,iBand)^2 + q_o{iBase,1}(iLayer,iUser,iBand)^2) / (b_o{iBase,1}(iLayer,iUser,iBand)) + ...
                                    (2 / b_o{iBase,1}(iLayer,iUser,iBand)) * (p_o{iBase,1}(iLayer,iUser,iBand) * (p(iLayer,iUser,iBand) - p_o{iBase,1}(iLayer,iUser,iBand))) + ...
                                    (2 / b_o{iBase,1}(iLayer,iUser,iBand)) * (q_o{iBase,1}(iLayer,iUser,iBand) * (q(iLayer,iUser,iBand) - q_o{iBase,1}(iLayer,iUser,iBand))) - ...
                                    (p_o{iBase,1}(iLayer,iUser,iBand)^2 + q_o{iBase,1}(iLayer,iUser,iBand)^2) / (2 * b_o{iBase,1}(iLayer,iUser,iBand)^2) * ...
                                    (b(iLayer,iUser,iBand) - b_o{iBase,1}(iLayer,iUser,iBand)) >= g(iLayer,iUser,iBand);
                                
                            end
                            
                        end
                        
                        norm(vec(M(:,:,:,iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                        
                        for iUser = 1:nUsers
                            for iLayer = 1:maxRank
                                intVector = [];
                                if SimStructs.userStruct{iUser,1}.baseNode ~= iBase
                                    for jUser = 1:kUsers
                                        for jLayer = 1:maxRank
                                            intVector = [intVector ; vW{iUser,1}(:,iLayer)' * cH{iBase,iBand}(:,:,iUser) * M(:,jLayer,jUser,iBand)];
                                        end
                                    end
                                    D{iLayer,iUser,iBand} : norm(intVector,2) <= fixedIF(iLayer,iUser,iBase,iBand);
                                end
                            end
                        end
                    end
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        norm(vec(t(:,iUser,:)),1) <= QueuedPkts(cUser,1) * log(2);
                    end
                    
                    cvx_end
                    
                    if strfind(cvx_status,'Solved')
                        cM{iBase,1} = M;
                        b_o{iBase,1} = b;
                        
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(:,iUser,:))) * log2(exp(1)),0);
                            SimParams.Debug.tempResource{2,1}{cUser,1} = [SimParams.Debug.tempResource{2,1}{cUser,1} sum(vec(t(:,iUser,:)))];
                            SimParams.Debug.tempResource{3,1}{cUser,1} = [SimParams.Debug.tempResource{3,1}{cUser,1} qDeviation];
                            for iBand = 1:nBands
                                SimParams.Debug.tempResource{4,1}{cUser,iBand} = [SimParams.Debug.tempResource{4,1}{cUser,iBand} sum(t(:,iUser,iBand))];
                            end
                        end
                        
                    else
                        b_o{iBase,1} = b_o{iBase,1} * 2;
                    end
                    
                    cD{iBase,1} = D;
                    
                end
                
                fixedIFH = fixedIF;
                fixedIF = zeros(maxRank,nUsers,nBases,nBands);
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:nUsers
                            baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                            if baseNode ~= iBase
                                for iLayer = 1:maxRank
                                    fixedIF(iLayer,iUser,iBase,iBand) = fixedIFH(iLayer,iUser,iBase,iBand) + alpha * (cD{baseNode,1}{iLayer,iUser,iBand} - cD{iBase,1}{iLayer,iUser,iBand});
                                end
                            end
                        end
                    end
                end
                
                display(fixedIF);
                totalCVXval = 0;
                for iBase = 1:nBases
                    for iUser = 1:nUsers
                        for iLayer = 1:maxRank
                            for iBand = 1:nBands
                                totalCVXval = totalCVXval + max(cD{iBase,1}{iLayer,iUser,iBand},0);
                            end
                        end
                    end
                end
                
                fixedIF = max(fixedIF,1e-2);
                if abs(totalCVXvalH - totalCVXval) <= 1e-2
                    masterIterate = 0;
                else
                    totalCVXvalH = totalCVXval;
                end
                
            end
            
            for iBand = 1:nBands
                for iBase = 1:nBases
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        currentH = cH{iBase,iBand}(:,:,cUser);
                        for iLayer = 1:maxRank
                            p_o{iBase,1}(iLayer,iUser,iBand) = real(vW{cUser,1}(:,iLayer)' * currentH * cM{iBase,1}(:,iLayer,iUser,iBand));
                            q_o{iBase,1}(iLayer,iUser,iBand) = imag(vW{cUser,1}(:,iLayer)' * currentH * cM{iBase,1}(:,iLayer,iUser,iBand));
                        end
                    end
                end
                
                for iBase = 1:nBases
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iLayer = 1:maxRank
                            R = SimParams.N * eye(SimParams.nRxAntenna);
                            for jBase = 1:nBases
                                for jUser = 1:usersPerCell(jBase,1)
                                    H = cH{jBase,iBand}(:,:,cUser);
                                    R = R + H * cM{jBase,1}(:,:,jUser,iBand) * cM{jBase,1}(:,:,jUser,iBand)' * H';
                                end
                            end
                            H = cH{iBase,iBand}(:,:,cUser);
                            vW{cUser,iBand}(:,iLayer) = R \ (H * cM{iBase,1}(:,iLayer,iUser,iBand));
                            vW{cUser,iBand}(:,iLayer) = vW{cUser,iBand}(:,iLayer) / norm(vW{cUser,iBand}(:,iLayer),2);
                        end
                    end
                end
            end
            
            if min(abs(cvx_optval - cvx_hist)) <= 1e-2
                reIterate = 0;
            else
                xIndex = xIndex + 1;
                cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
            end
            
        end
        
        updatePrecoders = 'false';
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = cM{iBase,1}(:,:,:,iBand);
            end
        end
        
end

if strcmp(updatePrecoders,'true')
    for iBase = 1:nBases
        for iBand = 1:nBands
            SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
            for iUser = 1:usersPerCell(iBase,1)
                cUser = cellUserIndices{iBase,1}(iUser,1);
                SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,iUser) = M(:,cUser,iBand);
            end
        end
    end
end

for iUser = 1:nUsers
    SimParams.Debug.tempResource{2,1}{iUser,1} = SimParams.Debug.tempResource{2,1}{iUser,1} * log2(exp(1));
    for iBand = 1:nBands
        SimParams.Debug.tempResource{4,1}{iUser,iBand} = SimParams.Debug.tempResource{4,1}{iUser,iBand} * log2(exp(1));
    end
end
