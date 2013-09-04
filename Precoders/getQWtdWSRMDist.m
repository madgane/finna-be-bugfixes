
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
    
    case 'DistJointAlloc'
       
        while masterIterate
            
            for iBase = 1:nBases
                
                xIndex = 0;
                cvx_hist = -500;
                localIterate = 1;
                p_o = ones(nUsers,nBands) / nUsers;
                q_o = ones(nUsers,nBands) / nUsers;
                b_o = ones(nUsers,nBands) * 10 + rand(nUsers,nBands);
                
                while localIterate
                    
                    cvx_begin
                    
                    dual varaibles dTO{nUsers,nBands} dFR{nUsers,nBands}
                    
                    expressions p(nUsers,nBands) q(nUsers,nBands)
                    variable M(SimParams.nTxAntenna,nUsers,nBands) complex
                    variables t(nUsers,nBands) b(nUsers,nBands) g(nUsers,nBands)
                    variables userObjective(nUsers,1) epiObjective
                    
                    minimize(epiObjective)
                    
                    subject to
                    
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        abs(QueuedPkts(cUser,1) - sum(t(cUser,:))) <= userObjective(cUser,1);
                    end
                                        
                    epiObjective >= norm(userObjective,1);
                    
                    for iBand = 1:nBands
                        
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            
                            intVector = sqrt(SimParams.N);
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            for jBase = 1:nBases
                                if jBase ~= iBase
                                    intVector = [intVector ; sqrt(fixedIF(cUser,jBase,iBand))];
                                end
                            end                            
                            
                            for jUser = 1:usersPerCell(iBase,1)
                                rUser = cellUserIndices{iBase,1}(iUser,1);
                                if rUser ~= cUser
                                    intVector = [intVector ; currentH * M(:,rUser,iBand)];
                                end
                            end
                            
                            dFR{cUser,iBand} : norm(intVector,2) <= sqrt(b(cUser,iBand));
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
                        
                        for jUser = 1:nUsers
                            if isempty(find(jUser == cellUserIndices{iBase,1}))
                                
                                intVector = [];
                                interH = cH{iBase,iBand}(:,:,jUser);
                                
                                for iUser = 1:usersPerCell(iBase,1)
                                    cUser = cellUserIndices{iBase,1}(iUser,1);
                                    intVector = [intVector ; interH * M(:,cUser,iBand)];
                                end
                                
                                dTO{jUser,iBand} : norm(intVector,2) <= sqrt(fixedIF(jUser,iBase,iBand));
                            end
                        end
                    end
                    
                    norm(vec(M(:,cellUserIndices{iBase,1},iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        norm(t(cUser,:),1) <= QueuedPkts(cUser,1) * log(2);
                    end
                    
                    cvx_end
                    
                    if strfind(cvx_status,'Solved')
                        
                        b_o = b;
                        M = full(M);
                        for iBand = 1:nBands
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
                        
                        if min(abs(cvx_optval - cvx_hist)) <= 1e-4
                            localIterate = 0;
                        else
                            xIndex = xIndex + 1;
                            cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                        end
                    else
                        b_o = b_o * 2;
                    end
                    
                end                
                
                cM{iBase,1} = M;
                cTO{iBase,1} = dTO;
                cFR{iBase,1} = dFR;
                                
            end
            
            fixedIFH = fixedIF;
            fixedIF = zeros(nUsers,nBases,nBands);
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                    for iBase = 1:nBases
                        if baseNode ~= iBase
                            fixedIF(iUser,iBase,iBand) = fixedIFH(iUser,iBase,iBand) + (cFR{baseNode,1}{iUser,iBand} - cTO{iBase,1}{iUser,iBand});
                        end
                    end                    
                end                
            end
            
            fixedIF = max(fixedIF,0);
            if norm(vec(fixedIF - fixedIFH),1) <= 1e-1
                masterIterate = 0;
            end
            
            norm(vec(fixedIF - fixedIFH),1)
            
        end
        
        updatePrecoders = 'false';
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,iUser) = cM{iBase,1}(:,cUser,iBand);
                end
            end
        end

                
    case 'DistGenAlloc'
        
        masterConvergence = 1;
        vW = cell(nUsers,nBands);
        maxRank = SimParams.maxRank;
        rankUsers = nUsers * maxRank;
        
        p_o = cell(nBases,1);
        q_o = cell(nBases,1);
        b_o = cell(nBases,1);
        lambda = ones(maxRank,nUsers,nBands);
        
        g_b = zeros(maxRank,nUsers,nBands);
        g_inter_cell = zeros(maxRank,nUsers,nBands);
        g_intra_cell = zeros(maxRank,nUsers,nBands);
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                vW{iBase,iBand} = zeros(SimParams.nTxAntenna,maxRank,usersPerCell(iBase,1));
            end
        end
        
        while masterConvergence
            
            xIndex = 0;
            reIterate = 1;
            cvx_hist = -500 * ones(2,1);
            firstAttempt = ones(nBases,1);
            
            for iBase = 1:nBases
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        vW{cUser,iBand} = ones(SimParams.nRxAntenna,maxRank) / sqrt(SimParams.nRxAntenna);
                    end
                end
            end
            
            for iBase = 1:nBases
                
                kUsers = usersPerCell(iBase,1);
                
                if firstAttempt(iBase,1)
                    firstAttempt(iBase,1) = 0;
                    p_o{iBase,1} = ones(maxRank,kUsers,nBands) / rankUsers;
                    q_o{iBase,1} = ones(maxRank,kUsers,nBands) / rankUsers;
                    b_o{iBase,1} = ones(maxRank,kUsers,nBands) * 10 + rand(maxRank,kUsers,nBands);
                else
                    for iBand = 1:nBands
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                R = SimParams.N * eye(SimParams.nRxAntenna);
                                for jBase = 1:nBases
                                    for jUser = 1:usersPerCell(jBase,1)
                                        rUser = cellUserIndices{jBase,1}(jUser,1);
                                        H = cH{jBase,iBand}(:,:,cUser);
                                        R = R + H * vM{rUser,iBand} * vM{rUser,iBand}' * H';
                                    end
                                end
                                H = cH{iBase,iBand}(:,:,cUser);
                                vW{cUser,iBand}(:,iLayer) = R \ (H * vM{cUser,iBand}(:,iLayer));
                                vW{cUser,iBand}(:,iLayer) = vW{cUser,iBand}(:,iLayer) / norm(vW{cUser,iBand}(:,iLayer),2);
                            end
                        end
                    end
                end
                
                while reIterate
                    
                    cvx_begin
                    
                    expressions p(maxRank,kUsers,nBands) q(maxRank,kUsers,nBands) expvar
                    variable M(SimParams.nTxAntenna,maxRank,kUsers,nBands) complex
                    variables t(maxRank,kUsers,nBands) b(maxRank,nUsers,nBands) g(maxRank,kUsers,nBands) inter_cell(maxRank,nUsers,nBands) intra_cell(maxRank,nUsers,nBands) interf(nBands,1)
                    variables userObjective(kUsers,1) epiObjective
                    
                    minimize(epiObjective)
                    
                    subject to
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        abs(QueuedPkts(cUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1);
                    end
                    
                    epiObjective >= norm(userObjective,1) + norm(interf,1);
                    
                    for iBand = 1:nBands
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                intVector = sqrt(SimParams.N);
                                currentH = cH{iBase,iBand}(:,:,cUser);
                                for jUser = 1:usersPerCell(iBase,1)
                                    rUser = cellUserIndices{iBase,1}(jUser,1);
                                    if rUser ~= cUser
                                        for jLayer = 1:maxRank
                                            intVector = [intVector ; vW{cUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,jUser,iBand)];
                                        end
                                    else
                                        for jLayer = 1:maxRank
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; vW{cUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,jUser,iBand)];
                                            end
                                        end
                                    end
                                end
                                
                                norm(intVector,2) <= sqrt(inter_cell(iLayer,cUser,iBand));
                                log(1 + g(iLayer,iUser,iBand)) >= t(iLayer,iUser,iBand);
                                
                                currentH = cH{iBase,iBand}(:,:,cUser);
                                p(iLayer,iUser,iBand) = real(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,iUser,iBand));
                                q(iLayer,iUser,iBand) = imag(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,iUser,iBand));
                                
                                q(iLayer,iUser,iBand) == 0;
                                
                                (p_o{iBase,1}(iLayer,iUser,iBand)^2 + q_o{iBase,1}(iLayer,iUser,iBand)^2) / (b_o{iBase,1}(iLayer,iUser,iBand)) + ...
                                    (2 / b_o{iBase,1}(iLayer,iUser,iBand)) * (p_o{iBase,1}(iLayer,iUser,iBand) * (p(iLayer,iUser,iBand) - p_o{iBase,1}(iLayer,iUser,iBand))) + ...
                                    (2 / b_o{iBase,1}(iLayer,iUser,iBand)) * (q_o{iBase,1}(iLayer,iUser,iBand) * (q(iLayer,iUser,iBand) - q_o{iBase,1}(iLayer,iUser,iBand))) - ...
                                    (p_o{iBase,1}(iLayer,iUser,iBand)^2 + q_o{iBase,1}(iLayer,iUser,iBand)^2) / (2 * b_o{iBase,1}(iLayer,iUser,iBand)^2) * ...
                                    (b(iLayer,cUser,iBand) - b_o{iBase,1}(iLayer,iUser,iBand)) >= g(iLayer,iUser,iBand);
                                
                            end
                            
                        end
                        
                        norm(vec(M(:,:,:,iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                        
                        for jBase = 1:nBases
                            if jBase ~= iBase
                                for jUser = 1:usersPerCell(jBase,1)
                                    cUser = cellUserIndices{jBase,1}(jUser,1);
                                    cHannel = cH{iBase,iBand}(:,:,cUser);
                                    for iLayer = 1:maxRank
                                        intVector = [];
                                        for kUser = 1:kUsers
                                            intVector = [intVector ; vec(vW{cUser,1}(:,iLayer)' * cHannel * M(:,:,kUser))];
                                        end
                                        sqrt(intra_cell(iLayer,cUser,iBand)) >= norm(intVector,2);
                                    end
                                end
                            end
                        end
                        
                        expvar = 0;
                        for kUser = 1:nUsers
                            if isempty(find(kUser == cellUserIndices{iBase,1}))
                                for iLayer = 1:maxRank
                                    expvar = expvar + intra_cell(iLayer,kUser,iBand) * lambda(iLayer,kUser,iBand);
                                end
                            else
                                for iLayer = 1:maxRank
                                    expvar = expvar + (inter_cell(iLayer,kUser,iBand) - b(iLayer,kUser,iBand)) * lambda(iLayer,kUser,iBand);
                                end
                            end
                        end
                        
                        interf(iBand,1) >= expvar;
                        
                    end
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        norm(vec(t(:,iUser,:)),1) <= QueuedPkts(cUser,1) * log(2);
                    end
                    
                    cvx_end
                    
                    if strfind(cvx_status,'Solved')
                        
                        for iBand = 1:nBands
                            for iUser = 1:usersPerCell(iBase,1)
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                currentH = cH{iBase,iBand}(:,:,cUser);
                                for iLayer = 1:maxRank
                                    b_o{iBase,1}(iLayer,iUser,iBand) = b(iLayer,cUser,iBand);
                                    p_o{iBase,1}(iLayer,iUser,iBand) = real(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,iUser,iBand));
                                    q_o{iBase,1}(iLayer,iUser,iBand) = imag(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,iUser,iBand));
                                    
                                    g_b(iLayer,cUser,iBand) = b(iLayer,cUser,iBand);
                                    g_inter_cell(iLayer,cUser,iBand) = inter_cell(iLayer,cUser,iBand);
                                    g_intra_cell(iLayer,cUser,iBand) = intra_cell(iLayer,cUser,iBand);
                                end
                                
                                vM{cUser,iBand} = M(:,:,iUser,iBand);
                                
                                if iBand == 1
                                    qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(:,cUser,:))) * log2(exp(1)),0);
                                    SimParams.Debug.tempResource{2,1}{cUser,1} = [SimParams.Debug.tempResource{2,1}{cUser,1} sum(vec(t(:,cUser,:)))];
                                    SimParams.Debug.tempResource{3,1}{cUser,1} = [SimParams.Debug.tempResource{3,1}{cUser,1} qDeviation];
                                end
                                SimParams.Debug.tempResource{4,1}{cUser,iBand} = [SimParams.Debug.tempResource{4,1}{cUser,iBand} sum(vec(t(:,cUser,iBand)))];
                            end
                        end
                        
                        if min(abs(cvx_optval - cvx_hist)) <= 1e-4
                            reIterate = 0;
                        else
                            xIndex = xIndex + 1;
                            cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                        end
                    else
                        b_o{iBase,1} = b_o{iBase,1} * 2;
                    end
                    
                end
                
            end
            
            lambda_old = lambda;
            lambda = lambda_old + 0.5 * (g_intra_cell + g_inter_cell - g_b);
            
            if norm(vec(lambda_old - lambda),2) <= 1e-2
                masterConvergence = 0;
            end
            
        end
        
        updatePrecoders = 'false';
        for iBase = 1:nBases
            for iBand = 1:nBands
                P = [];
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P vM{iBase,iBand}(:,:,cUser)];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
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
