
function [SimParams,SimStructs] = getQWtdWSRM(SimParams,SimStructs)

cH = SimStructs.linkChan;
nBases = SimParams.nBases;
nBands = SimParams.nBands;

updatePrecoders = 'true';
usersPerCell = zeros(nBases,1);
cellUserIndices = cell(nBases,1);

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
    
    case 'BandAlloc'
        
        currentDesign = SimParams.PrecodingMethod;
        currentApproach = SimParams.weightedSumRateMethod;
        
        SimParams.PrecodingMethod = 'Best_WMMSE_Method';
        SimParams.weightedSumRateMethod = 'PerformScheduling';
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    SimStructs.userStruct{cUser,1}.weighingFactor = QueuedPkts(iUser,1);
                end
            end
            
            SimParams.Debug.privateExchanges.takeOverBand = iBand;
            [SimParams, SimStructs] = getWeightedMMSEDesign(SimParams,SimStructs);
            
            [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,iBand);
            QueuedPkts = max(QueuedPkts - SimParams.Debug.privateExchanges.resAllocation(iBand,:)',0);
            
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
        
    case 'ConstrBandAlloc'
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    SimStructs.userStruct{cUser,1}.weighingFactor = QueuedPkts(iUser,1);
                end
            end
            
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
                
                maximize(epiObjective)
                
                subject to
                
                for iUser = 1:nUsers
                    QueuedPkts(iUser,1) * t(iUser,1) >= userObjective(iUser,1);
                end
                
                epiObjective <= sum(userObjective,1);
                
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
            
        end
        
        updatePrecoders = 'false';
        
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

end
