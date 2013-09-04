
function [SimParams,SimStructs] = getQWtdWSRM(SimParams,SimStructs)

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
bandRateMax = zeros(SimParams.nUsers,nBands);

for iBase = 1:nBases
    for iUser = 1:usersPerCell(iBase,1)
        cUser = cellUserIndices{iBase,1}(iUser,1);
        QueuedPkts(cUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
    end
end

switch SimParams.weightedSumRateMethod
    
    case 'BandAlloc'
        
        currentDesign = SimParams.PrecodingMethod;
        currentApproach = SimParams.weightedSumRateMethod;
        
        SimParams.PrecodingMethod = 'Best_WMMSE_Method';
        SimParams.weightedSumRateMethod = 'PerformScheduling';
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    SimStructs.userStruct{cUser,1}.weighingFactor = QueuedPkts(cUser,1);
                end
            end
            
            SimParams.Debug.privateExchanges.takeOverBand = iBand;
            [SimParams, SimStructs] = getWeightedMMSEDesign(SimParams,SimStructs);
            
            [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,iBand);
            QueuedPkts = max(QueuedPkts - SimParams.Debug.privateExchanges.resAllocation(iBand,:)',0);
            
        end
        
        for iBase = 1:nBases
            for iUser = 1:usersPerCell(iBase,1)
                cUser = cellUserIndices{iBase,1}(iUser,1);
                sumRateOverBand = [];bandMaxRate = 0;
                for iBand = 1:nBands
                    sumRateOverBand = [sumRateOverBand (SimParams.Debug.tempResource{4,1}{cUser,iBand} + bandMaxRate)];
                    bandMaxRate = max(sumRateOverBand);
                end
                SimParams.Debug.tempResource{2,1}{cUser,1} = sumRateOverBand;
            end
        end
        
        
        SimParams.PrecodingMethod = currentDesign;
        SimParams.weightedSumRateMethod = currentApproach;
        SimParams.privateExchanges = rmfield(SimParams.Debug.privateExchanges,'takeOverBand');
        
        updatePrecoders = 'false';
        
    case 'JointAlloc'
        
        reIterate = 1;
        cvx_hist = -500 * ones(2,1);
        
        xIndex = 0;
        p_o = ones(nUsers,nBands) / nUsers;
        q_o = ones(nUsers,nBands) / nUsers;
        b_o = ones(nUsers,nBands) * 10 + rand(nUsers,nBands);
        
        while reIterate
            
            cvx_begin
            
            expressions p(nUsers,nBands) q(nUsers,nBands)
            variable M(SimParams.nTxAntenna,nUsers,nBands) complex
            variables t(nUsers,nBands) b(nUsers,nBands) g(nUsers,nBands)
            variables userObjective(nUsers,1) epiObjective
            
            minimize(epiObjective)
            
            subject to
            
            for iUser = 1:nUsers
                abs(QueuedPkts(iUser,1) - sum(t(iUser,:))) <= userObjective(iUser,1);
            end
            
            epiObjective >= norm(userObjective,1);
            
            for iBase = 1:nBases
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        intVector = sqrt(SimParams.N);
                        
                        for jBase = 1:nBases
                            currentH = cH{jBase,iBand}(:,:,cUser);
                            for jUser = 1:usersPerCell(jBase,1)
                                rUser = cellUserIndices{jBase,1}(jUser,1);
                                if rUser ~= cUser
                                    intVector = [intVector ; currentH * M(:,rUser,iBand)];
                                end
                            end
                        end
                        
                        norm(intVector,2) <= sqrt(b(cUser,iBand));
                        log(1 + g(cUser,iBand)) >= t(cUser,iBand);
                        
                        currentH = cH{iBase,iBand}(:,:,cUser);
                        p(cUser,iBand) = real(currentH * M(:,cUser,iBand));
                        q(cUser,iBand) = imag(currentH * M(:,cUser,iBand));
                        
                        q(cUser,iBand) == 0;
                        
                        (p_o(cUser,iBand)^2 + q_o(cUser,iBand)^2) / (b_o(cUser,iBand)) + ...
                            (2 / b_o(cUser,iBand)) * (p_o(cUser,iBand) * (p(cUser,iBand) - p_o(cUser,iBand))) + ...
                            (2 / b_o(cUser,iBand)) * (q_o(cUser,iBand) * (q(cUser,iBand) - q_o(cUser,iBand))) - ...
                            (p_o(cUser,iBand)^2 + q_o(cUser,iBand)^2) / (2 * b_o(cUser,iBand)^2) * ...
                            (b(cUser,iBand) - b_o(cUser,iBand)) >= g(cUser,iBand);
                        
                    end
                    
                    norm(vec(M(:,cellUserIndices{iBase,1},iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    
                end
                
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    norm(t(cUser,:),1) <= QueuedPkts(cUser,1) * log(2);
                end
                
            end
            
            cvx_end
            
            if strfind(cvx_status,'Solved')
                
                b_o = b;
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            p_o(cUser,iBand) = real(currentH * M(:,cUser,iBand));
                            q_o(cUser,iBand) = imag(currentH * M(:,cUser,iBand));
                            
                            if iBand == 1
                                qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(cUser,:))) * log2(exp(1)),0);
                                SimParams.Debug.tempResource{2,1}{cUser,1} = [SimParams.Debug.tempResource{2,1}{cUser,1} sum(vec(t(cUser,:)))];
                                SimParams.Debug.tempResource{3,1}{cUser,1} = [SimParams.Debug.tempResource{3,1}{cUser,1} qDeviation];
                            end
                            SimParams.Debug.tempResource{4,1}{cUser,iBand} = [SimParams.Debug.tempResource{4,1}{cUser,iBand} t(cUser,iBand)];
                        end
                    end
                end
                
                if min(abs(cvx_optval - cvx_hist)) <= 1e-4
                    reIterate = 0;
                else
                    xIndex = xIndex + 1;
                    cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                end
            else
                b_o = b_o * 2;
            end
            
        end
        
    case 'GMApproxAlloc'
        
        reIterate = 1;
        cvx_hist = -500 * ones(2,1);
        
        xIndex = 0;
        p_o = ones(nUsers,nBands) / nUsers;
        q_o = ones(nUsers,nBands) / nUsers;
        phi = ones(nUsers,nBands) * 10 + rand(nUsers,nBands);
        
        while reIterate
            
            cvx_begin
            
            expressions p(nUsers,nBands) q(nUsers,nBands)
            variable M(SimParams.nTxAntenna,nUsers,nBands) complex
            variables t(nUsers,nBands) b(nUsers,nBands) g(nUsers,nBands)
            variables userObjective(nUsers,1) epiObjective
            
            minimize(epiObjective)
            
            subject to
            
            for iUser = 1:nUsers
                abs(QueuedPkts(iUser,1) - sum(t(iUser,:))) <= userObjective(iUser,1);
            end
            
            epiObjective >= norm(userObjective,1);
            
            for iBase = 1:nBases
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        intVector = sqrt(SimParams.N);
                        
                        for jBase = 1:nBases
                            currentH = cH{jBase,iBand}(:,:,cUser);
                            for jUser = 1:usersPerCell(jBase,1)
                                rUser = cellUserIndices{jBase,1}(jUser,1);
                                if rUser ~= cUser
                                    intVector = [intVector ; currentH * M(:,rUser,iBand)];
                                end
                            end
                        end
                        
                        norm(intVector,2) <= sqrt(b(cUser,iBand));
                        log(g(cUser,iBand)) >= t(cUser,iBand);
                        
                        currentH = cH{iBase,iBand}(:,:,cUser);
                        p(cUser,iBand) = real(currentH * M(:,cUser,iBand));
                        q(cUser,iBand) = imag(currentH * M(:,cUser,iBand));
                        
                        q(cUser,iBand) == 0;
                        
                        p_o(cUser,iBand)^2 + q_o(cUser,iBand)^2 + ...
                            2 * (p_o(cUser,iBand) * (p(cUser,iBand) - p_o(cUser,iBand))) + ...
                            2 * (q_o(cUser,iBand) * (q(cUser,iBand) - q_o(cUser,iBand))) + b(cUser,iBand) >= ...
                            0.5 * (phi(cUser,iBand) * g(cUser,iBand)^2 + (1 / phi(cUser,iBand)) * b(cUser,iBand)^2);
                        
                    end
                    
                    norm(vec(M(:,cellUserIndices{iBase,1},iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    
                end
                
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    norm(t(cUser,:),1) <= QueuedPkts(cUser,1) * log(2);
                end
                
            end
            
            cvx_end
            
            if strfind(cvx_status,'Solved')
                
                phi = b ./ g;
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            p_o(cUser,iBand) = real(currentH * M(:,cUser,iBand));
                            q_o(cUser,iBand) = imag(currentH * M(:,cUser,iBand));
                            
                            if iBand == 1
                                qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(cUser,:))) * log2(exp(1)),0);
                                SimParams.Debug.tempResource{2,1}{cUser,1} = [SimParams.Debug.tempResource{2,1}{cUser,1} sum(vec(t(cUser,:)))];
                                SimParams.Debug.tempResource{3,1}{cUser,1} = [SimParams.Debug.tempResource{3,1}{cUser,1} qDeviation];
                            end
                            
                            SimParams.Debug.tempResource{4,1}{cUser,iBand} = [SimParams.Debug.tempResource{4,1}{cUser,iBand} t(cUser,iBand)];
                            
                        end
                    end
                end
                
                if min(abs(cvx_optval - cvx_hist)) <= 1e-4
                    reIterate = 0;
                else
                    xIndex = xIndex + 1;
                    cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                end
            else
                phi = ones(nUsers,nBands) * 10 + rand(nUsers,nBands);
            end
            
        end
        
    case 'JointBandAlloc'
        
        for iBand = 1:nBands
            
            reIterate = 1;
            cvx_hist = -500 * ones(2,1);
            
            xIndex = 0;
            p_o = ones(nUsers,1) / nUsers;
            q_o = ones(nUsers,1) / nUsers;
            b_o = ones(nUsers,1) * 10 + rand(nUsers,1);
            
            while reIterate
                
                cvx_begin
                
                expressions p(nUsers,1) q(nUsers,1)
                variable M(SimParams.nTxAntenna,nUsers) complex
                variables t(nUsers,1) b(nUsers,1) g(nUsers,1)
                variables userObjective(nUsers,1) epiObjective
                
                minimize(epiObjective)
                
                subject to
                
                for iUser = 1:nUsers
                    abs(QueuedPkts(iUser,1) - t(iUser,1)) <= userObjective(iUser,1);
                end
                
                epiObjective >= norm(userObjective,1);
                
                for iBase = 1:nBases
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        intVector = sqrt(SimParams.N);
                        
                        for jBase = 1:nBases
                            currentH = cH{jBase,iBand}(:,:,cUser);
                            for jUser = 1:usersPerCell(jBase,1)
                                rUser = cellUserIndices{jBase,1}(jUser,1);
                                if rUser ~= cUser
                                    intVector = [intVector ; currentH * M(:,rUser)];
                                end
                            end
                        end
                        
                        norm(intVector,2) <= sqrt(b(cUser,1));
                        log(1 + g(cUser,1)) >= t(cUser,1);
                        
                        currentH = cH{iBase,iBand}(:,:,cUser);
                        p(cUser,1) = real(currentH * M(:,cUser));
                        q(cUser,1) = imag(currentH * M(:,cUser));
                        
                        q(cUser,1) == 0;
                        
                        (p_o(cUser,1)^2 + q_o(cUser,1)^2) / (b_o(cUser,1)) + ...
                            (2 / b_o(cUser,1)) * (p_o(cUser,1) * (p(cUser,1) - p_o(cUser,1))) + ...
                            (2 / b_o(cUser,1)) * (q_o(cUser,1) * (q(cUser,1) - q_o(cUser,1))) - ...
                            (p_o(cUser,1)^2 + q_o(cUser,1)^2) / (2 * b_o(cUser,1)^2) * ...
                            (b(cUser,1) - b_o(cUser,1)) >= g(cUser,1);
                        
                    end
                    
                    norm(vec(M(:,cellUserIndices{iBase,1})),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        norm(t(cUser,1),1) <= QueuedPkts(cUser,1) * log(2);
                    end
                    
                end
                
                cvx_end
                
                if strfind(cvx_status,'Solved')
                    
                    b_o = b;
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            p_o(cUser,1) = real(currentH * M(:,cUser));
                            q_o(cUser,1) = imag(currentH * M(:,cUser));
                            
                            qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(cUser,:))) * log2(exp(1)),0);
                            SimParams.Debug.tempResource{2,1}{cUser,1} = [SimParams.Debug.tempResource{2,1}{cUser,1} sum(vec(t(cUser,:))) + sum(bandRateMax(cUser,1:(iBand-1)))];
                            SimParams.Debug.tempResource{3,1}{cUser,1} = [SimParams.Debug.tempResource{3,1}{cUser,1} qDeviation];
                            SimParams.Debug.tempResource{4,1}{cUser,iBand} = [SimParams.Debug.tempResource{4,1}{cUser,iBand} t(cUser,1)];
                        end
                    end
                    
                    if min(abs(cvx_optval - cvx_hist)) <= 1e-4
                        reIterate = 0;
                    else
                        xIndex = xIndex + 1;
                        cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                    end
                else
                    b_o = b_o * 2;
                end
                
            end
            
            for iBase = 1:nBases
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,iUser) = M(:,cUser);
                end
            end
            
            [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,iBand);
            QueuedPkts = max(QueuedPkts - SimParams.Debug.privateExchanges.resAllocation(iBand,:)',0);
            bandRateMax(:,iBand) = log(2) * SimParams.Debug.privateExchanges.resAllocation(iBand,:)';
        end
        
        updatePrecoders = 'false';
        
    case 'GenAlloc'
        
        vW = cell(nUsers,nBands);
        maxRank = SimParams.maxRank;
        rankUsers = nUsers * maxRank;
        
        xIndex = 0;
        reIterate = 1;
        cvx_hist = -500 * ones(2,1);
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    vW{cUser,iBand} = ones(SimParams.nRxAntenna,maxRank) / sqrt(SimParams.nRxAntenna);
                end
            end
        end
        
        p_o = ones(maxRank,nUsers,nBands) / rankUsers;
        q_o = ones(maxRank,nUsers,nBands) / rankUsers;
        b_o = ones(maxRank,nUsers,nBands) * 10 + rand(maxRank,nUsers,nBands);
        
        while reIterate
            
            cvx_begin
            
            expressions p(maxRank,nUsers,nBands) q(maxRank,nUsers,nBands)
            variable M(SimParams.nTxAntenna,maxRank,nUsers,nBands) complex
            variables t(maxRank,nUsers,nBands) b(maxRank,nUsers,nBands) g(maxRank,nUsers,nBands)
            variables userObjective(nUsers,1) epiObjective
            
            minimize(epiObjective)
            
            subject to
            
            for iUser = 1:nUsers
                abs(QueuedPkts(iUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1);
            end
            
            epiObjective >= norm(userObjective,1);
            
            for iBase = 1:nBases
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iLayer = 1:maxRank
                            intVector = sqrt(SimParams.N);
                            
                            for jBase = 1:nBases
                                currentH = cH{jBase,iBand}(:,:,cUser);
                                for jUser = 1:usersPerCell(jBase,1)
                                    rUser = cellUserIndices{jBase,1}(jUser,1);
                                    if rUser ~= cUser
                                        for jLayer = 1:maxRank
                                            intVector = [intVector ; vW{cUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,rUser,iBand)];
                                        end
                                    else
                                        for jLayer = 1:maxRank
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; vW{cUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,rUser,iBand)];
                                            end
                                        end
                                    end
                                end
                            end
                            
                            norm(intVector,2) <= sqrt(b(iLayer,cUser,iBand));
                            log(1 + g(iLayer,cUser,iBand)) >= t(iLayer,cUser,iBand);
                            
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            p(iLayer,cUser,iBand) = real(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                            q(iLayer,cUser,iBand) = imag(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                            
                            q(iLayer,cUser,iBand) == 0;
                            
                            (p_o(iLayer,cUser,iBand)^2 + q_o(iLayer,cUser,iBand)^2) / (b_o(iLayer,cUser,iBand)) + ...
                                (2 / b_o(iLayer,cUser,iBand)) * (p_o(iLayer,cUser,iBand) * (p(iLayer,cUser,iBand) - p_o(iLayer,cUser,iBand))) + ...
                                (2 / b_o(iLayer,cUser,iBand)) * (q_o(iLayer,cUser,iBand) * (q(iLayer,cUser,iBand) - q_o(iLayer,cUser,iBand))) - ...
                                (p_o(iLayer,cUser,iBand)^2 + q_o(iLayer,cUser,iBand)^2) / (2 * b_o(iLayer,cUser,iBand)^2) * ...
                                (b(iLayer,cUser,iBand) - b_o(iLayer,cUser,iBand)) >= g(iLayer,cUser,iBand);
                            
                        end
                        
                    end
                    
                    norm(vec(M(:,:,cellUserIndices{iBase,1},iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    
                end
                
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    norm(vec(t(:,cUser,:)),1) <= QueuedPkts(cUser,1) * log(2);
                end
                
            end
            
            cvx_end
            
            if strfind(cvx_status,'Solved')
                
                b_o = b;
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            for iLayer = 1:maxRank
                                p_o(iLayer,cUser,iBand) = real(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                                q_o(iLayer,cUser,iBand) = imag(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                            end
                            
                            if iBand == 1
                                qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(:,cUser,:))) * log2(exp(1)),0);
                                SimParams.Debug.tempResource{2,1}{cUser,1} = [SimParams.Debug.tempResource{2,1}{cUser,1} sum(vec(t(:,cUser,:)))];
                                SimParams.Debug.tempResource{3,1}{cUser,1} = [SimParams.Debug.tempResource{3,1}{cUser,1} qDeviation];
                            end
                            
                            SimParams.Debug.tempResource{4,1}{cUser,iBand} = [SimParams.Debug.tempResource{4,1}{cUser,iBand} sum(vec(t(:,cUser,iBand)))];
                            
                        end
                    end
                end
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                R = SimParams.N * eye(SimParams.nRxAntenna);
                                for jBase = 1:nBases
                                    for jUser = 1:usersPerCell(jBase,1)
                                        rUser = cellUserIndices{jBase,1}(jUser,1);
                                        H = cH{jBase,iBand}(:,:,cUser);
                                        R = R + H * M(:,:,rUser,iBand) * M(:,:,rUser,iBand)' * H';
                                    end
                                end
                                H = cH{iBase,iBand}(:,:,cUser);
                                vW{cUser,iBand}(:,iLayer) = R \ (H * M(:,iLayer,cUser,iBand));
                                vW{cUser,iBand}(:,iLayer) = vW{cUser,iBand}(:,iLayer) / norm(vW{cUser,iBand}(:,iLayer),2);
                            end
                        end
                    end
                end
                
                if min(abs(cvx_optval - cvx_hist)) <= 1e-4
                    reIterate = 0;
                else
                    xIndex = xIndex + 1;
                    cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                end
            else
                b_o = b_o * 2;
            end
            
        end
        
        updatePrecoders = 'false';
        for iBase = 1:nBases
            for iBand = 1:nBands
                P = [];
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P M(:,:,cUser,iBand)];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
        end
        
    case 'GenBandAlloc'
        
        for iBand = 1:nBands
            
            vW = cell(nUsers,1);
            maxRank = SimParams.maxRank;
            rankUsers = nUsers * maxRank;
            
            xIndex = 0;
            reIterate = 1;
            cvx_hist = -500 * ones(2,1);
            
            for iBase = 1:nBases
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    vW{cUser,1} = ones(SimParams.nRxAntenna,maxRank) / sqrt(SimParams.nRxAntenna);
                end
            end
            
            
            p_o = ones(maxRank,nUsers) / rankUsers;
            q_o = ones(maxRank,nUsers) / rankUsers;
            b_o = ones(maxRank,nUsers) * 10 + rand(maxRank,nUsers);
            
            while reIterate
                
                cvx_begin
                
                expressions p(maxRank,nUsers) q(maxRank,nUsers)
                variable M(SimParams.nTxAntenna,maxRank,nUsers) complex
                variables t(maxRank,nUsers) b(maxRank,nUsers) g(maxRank,nUsers)
                variables userObjective(nUsers,1) epiObjective
                
                minimize(epiObjective)
                
                subject to
                
                for iUser = 1:nUsers
                    abs(QueuedPkts(iUser,1) - sum(vec(t(:,iUser)))) <= userObjective(iUser,1);
                end
                
                epiObjective >= norm(userObjective,1);
                
                for iBase = 1:nBases
                    for iUser = 1:usersPerCell(iBase,1)
                        
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iLayer = 1:maxRank
                            intVector = sqrt(SimParams.N);
                            
                            for jBase = 1:nBases
                                currentH = cH{jBase,iBand}(:,:,cUser);
                                for jUser = 1:usersPerCell(jBase,1)
                                    rUser = cellUserIndices{jBase,1}(jUser,1);
                                    if rUser ~= cUser
                                        for jLayer = 1:maxRank
                                            intVector = [intVector ; vW{cUser,1}(:,iLayer)' * currentH * M(:,jLayer,rUser)];
                                        end
                                    else
                                        for jLayer = 1:maxRank
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; vW{cUser,1}(:,iLayer)' * currentH * M(:,jLayer,rUser)];
                                            end
                                        end
                                    end
                                end
                            end
                            
                            norm(intVector,2) <= sqrt(b(iLayer,cUser));
                            log(1 + g(iLayer,cUser)) >= t(iLayer,cUser);
                            
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            p(iLayer,cUser) = real(vW{cUser,1}(:,iLayer)' * currentH * M(:,iLayer,cUser));
                            q(iLayer,cUser) = imag(vW{cUser,1}(:,iLayer)' * currentH * M(:,iLayer,cUser));
                            
                            q(iLayer,cUser) == 0;
                            
                            (p_o(iLayer,cUser)^2 + q_o(iLayer,cUser)^2) / (b_o(iLayer,cUser)) + ...
                                (2 / b_o(iLayer,cUser)) * (p_o(iLayer,cUser) * (p(iLayer,cUser) - p_o(iLayer,cUser))) + ...
                                (2 / b_o(iLayer,cUser)) * (q_o(iLayer,cUser) * (q(iLayer,cUser) - q_o(iLayer,cUser))) - ...
                                (p_o(iLayer,cUser)^2 + q_o(iLayer,cUser)^2) / (2 * b_o(iLayer,cUser)^2) * ...
                                (b(iLayer,cUser) - b_o(iLayer,cUser)) >= g(iLayer,cUser);
                            
                        end
                        
                    end
                    
                    norm(vec(M(:,:,cellUserIndices{iBase,1})),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        norm(vec(t(:,cUser)),1) <= QueuedPkts(cUser,1) * log(2);
                    end
                    
                end
                
                cvx_end
                
                if strfind(cvx_status,'Solved')
                    
                    b_o = b;
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            for iLayer = 1:maxRank
                                p_o(iLayer,cUser) = real(vW{cUser,1}(:,iLayer)' * currentH * M(:,iLayer,cUser));
                                q_o(iLayer,cUser) = imag(vW{cUser,1}(:,iLayer)' * currentH * M(:,iLayer,cUser));
                            end
                            
                            qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(:,cUser,:))) * log2(exp(1)),0);
                            SimParams.Debug.tempResource{2,1}{cUser,1} = [SimParams.Debug.tempResource{2,1}{cUser,1} sum(vec(t(:,cUser,:))) + sum(bandRateMax(cUser,1:(iBand - 1)))];
                            SimParams.Debug.tempResource{3,1}{cUser,1} = [SimParams.Debug.tempResource{3,1}{cUser,1} qDeviation];
                            SimParams.Debug.tempResource{4,1}{cUser,iBand} = [SimParams.Debug.tempResource{4,1}{cUser,iBand} sum(vec(t(:,cUser,1)))];
                        end
                    end
                    
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                R = SimParams.N * eye(SimParams.nRxAntenna);
                                for jBase = 1:nBases
                                    for jUser = 1:usersPerCell(jBase,1)
                                        rUser = cellUserIndices{jBase,1}(jUser,1);
                                        H = cH{jBase,iBand}(:,:,cUser);
                                        R = R + H * M(:,:,rUser) * M(:,:,rUser)' * H';
                                    end
                                end
                                H = cH{iBase,iBand}(:,:,cUser);
                                vW{cUser,1}(:,iLayer) = R \ (H * M(:,iLayer,cUser));
                                vW{cUser,1}(:,iLayer) = vW{cUser,1}(:,iLayer) / norm(vW{cUser,1}(:,iLayer),2);
                            end
                        end
                    end
                    
                    if min(abs(cvx_optval - cvx_hist)) <= 1e-4
                        reIterate = 0;
                    else
                        xIndex = xIndex + 1;
                        cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                    end
                else
                    b_o = b_o * 2;
                end
                
            end
            
            updatePrecoders = 'false';
            for iBase = 1:nBases
                P = [];
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P M(:,:,cUser)];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
            
            [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,iBand);
            QueuedPkts = max(QueuedPkts - SimParams.Debug.privateExchanges.resAllocation(iBand,:)',0);
            bandRateMax(:,iBand) = log(2) * SimParams.Debug.privateExchanges.resAllocation(iBand,:)';
            
        end
        
    case 'GlobalAlloc'
        
        vW = cell(nUsers,nBands);
        maxRank = SimParams.maxRank;
        rankUsers = nUsers * maxRank;
        
        xIndex = 0;
        reIterate = 1;
        cvx_hist = -500 * ones(2,1);
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    vW{cUser,iBand} = ones(SimParams.nRxAntenna,maxRank) / sqrt(SimParams.nRxAntenna);
                end
            end
        end
        
        p_o = ones(maxRank,nUsers,nBands) / rankUsers;
        q_o = ones(maxRank,nUsers,nBands) / rankUsers;
        b_o = ones(maxRank,nUsers,nBands) * 10 + rand(maxRank,nUsers,nBands);
        
        while reIterate
            
            cvx_begin
            
            expressions p(maxRank,nUsers,nBands) q(maxRank,nUsers,nBands)
            variable M(SimParams.nTxAntenna,maxRank,nUsers,nBands) complex
            variables t(maxRank,nUsers,nBands) b(maxRank,nUsers,nBands) g(maxRank,nUsers,nBands)
            variables userObjective(nUsers,1) epiObjective
            
            minimize(epiObjective)
            
            subject to
            
            for iUser = 1:nUsers
                abs(QueuedPkts(iUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1);
            end
            
            epiObjective >= norm(userObjective,1);
            
            for iBase = 1:nBases
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iLayer = 1:maxRank
                            intVector = sqrt(SimParams.N);
                            
                            for jBase = 1:nBases
                                currentH = cH{jBase,iBand}(:,:,cUser);
                                for jUser = 1:usersPerCell(jBase,1)
                                    rUser = cellUserIndices{jBase,1}(jUser,1);
                                    if rUser ~= cUser
                                        for jLayer = 1:maxRank
                                            intVector = [intVector ; vW{cUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,rUser,iBand)];
                                        end
                                    else
                                        for jLayer = 1:maxRank
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; vW{cUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,rUser,iBand)];
                                            end
                                        end
                                    end
                                end
                            end
                            
                            norm(intVector,2) <= sqrt(b(iLayer,cUser,iBand));
                            log(1 + g(iLayer,cUser,iBand)) >= t(iLayer,cUser,iBand);
                            
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            p(iLayer,cUser,iBand) = real(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                            q(iLayer,cUser,iBand) = imag(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                            
                            q(iLayer,cUser,iBand) == 0;
                            
                            (p_o(iLayer,cUser,iBand)^2 + q_o(iLayer,cUser,iBand)^2) / (b_o(iLayer,cUser,iBand)) + ...
                                (2 / b_o(iLayer,cUser,iBand)) * (p_o(iLayer,cUser,iBand) * (p(iLayer,cUser,iBand) - p_o(iLayer,cUser,iBand))) + ...
                                (2 / b_o(iLayer,cUser,iBand)) * (q_o(iLayer,cUser,iBand) * (q(iLayer,cUser,iBand) - q_o(iLayer,cUser,iBand))) - ...
                                (p_o(iLayer,cUser,iBand)^2 + q_o(iLayer,cUser,iBand)^2) / (2 * b_o(iLayer,cUser,iBand)^2) * ...
                                (b(iLayer,cUser,iBand) - b_o(iLayer,cUser,iBand)) >= g(iLayer,cUser,iBand);
                            
                        end
                        
                    end
                    
                end
                
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    norm(vec(t(:,cUser,:)),1) <= QueuedPkts(cUser,1) * log(2);
                end
                
                norm(vec(M(:,:,cellUserIndices{iBase,1},:)),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower));
                
            end
            
            cvx_end
            
            if strfind(cvx_status,'Solved')
                
                b_o = b;
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            for iLayer = 1:maxRank
                                p_o(iLayer,cUser,iBand) = real(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                                q_o(iLayer,cUser,iBand) = imag(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                            end
                            
                            if iBand == 1
                                qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(:,cUser,:))) * log2(exp(1)),0);
                                SimParams.Debug.tempResource{2,1}{cUser,1} = [SimParams.Debug.tempResource{2,1}{cUser,1} sum(vec(t(:,cUser,:)))];
                                SimParams.Debug.tempResource{3,1}{cUser,1} = [SimParams.Debug.tempResource{3,1}{cUser,1} qDeviation];
                            end
                            SimParams.Debug.tempResource{4,1}{cUser,iBand} = [SimParams.Debug.tempResource{4,1}{cUser,iBand} sum(vec(t(:,cUser,iBand)))];
                        end
                    end
                end
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                R = SimParams.N * eye(SimParams.nRxAntenna);
                                for jBase = 1:nBases
                                    for jUser = 1:usersPerCell(jBase,1)
                                        rUser = cellUserIndices{jBase,1}(jUser,1);
                                        H = cH{jBase,iBand}(:,:,cUser);
                                        R = R + H * M(:,:,rUser,iBand) * M(:,:,rUser,iBand)' * H';
                                    end
                                end
                                H = cH{iBase,iBand}(:,:,cUser);
                                vW{cUser,iBand}(:,iLayer) = R \ (H * M(:,iLayer,cUser,iBand));
                                vW{cUser,iBand}(:,iLayer) = vW{cUser,iBand}(:,iLayer) / norm(vW{cUser,iBand}(:,iLayer),2);
                            end
                        end
                    end
                end
                
                if min(abs(cvx_optval - cvx_hist)) <= 1e-4
                    reIterate = 0;
                else
                    xIndex = xIndex + 1;
                    cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                end
            else
                b_o = b_o * 2;
            end
            
        end
        
        updatePrecoders = 'false';
        for iBase = 1:nBases
            for iBand = 1:nBands
                P = [];
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P M(:,:,cUser,iBand)];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
        end
        
        
    case 'PrimalJAlloc'
        
        alpha = 0.9;
        reIterate = 1;
        cvx_hist = -500 * ones(2,1);
        p_o = cell(nBases,1);q_o = cell(nBases,1);b_o = cell(nBases,1);cM = cell(nBases,1);
                
        xIndex = 0;
        for iBase = 1:nBases
            p_o{iBase,1} = ones(usersPerCell(iBase,1),nBands) / usersPerCell(iBase,1);
            q_o{iBase,1} = ones(usersPerCell(iBase,1),nBands) / usersPerCell(iBase,1);
            b_o{iBase,1} = ones(usersPerCell(iBase,1),nBands) * 10 + rand(usersPerCell(iBase,1),nBands);
        end
        
        while reIterate
        
            masterIterate = 1;totalCVXvalH = -50;
            fixedIF = ones(nUsers,nBases,nBands) * 10;
        
            while masterIterate
                
                totalCVXval = 0;
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
                                    intVector = [intVector ; sqrt(fixedIF(cUser,jBase,iBand))];
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
                                D{iUser,iBand} : norm(intVector,2) <= sqrt(fixedIF(iUser,iBase,iBand));
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
                    totalCVXval = totalCVXval + cvx_optval;                    
                    
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
                display(abs(totalCVXvalH - totalCVXval));
                
                fixedIF = max(fixedIF,0);
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
                        p_o{iBase,1}(iUser,iBand) = real(currentH * cM{iBase,1}(:,iUser,iBand));
                        q_o{iBase,1}(iUser,iBand) = imag(currentH * cM{iBase,1}(:,iUser,iBand));
                    end
                end
            end
            
            if min(abs(cvx_optval - cvx_hist)) <= 1e-3
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

end
