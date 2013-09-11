
function [SimParams,SimStructs] = getQWtdWSRMDist(SimParams,SimStructs)

cH = SimStructs.linkChan;
nBases = SimParams.nBases;
nBands = SimParams.nBands;

usersPerCell = zeros(nBases,1);
cellUserIndices = cell(nBases,1);
cellNeighbourIndices = cell(nBases,1);

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

for iBase = 1:nBases
    for jBase = 1:nBases
        if jBase ~= iBase
            cellNeighbourIndices{iBase,1} = [cellNeighbourIndices{iBase,1} ; cellUserIndices{jBase,1}];
        end
    end
end

switch SimParams.weightedSumRateMethod
    
    case 'PrimalMethod'
        
        alpha = 0.5;
        nLayers = SimParams.maxRank;
        mIterationsSCA = 50;mIterationsSG = 20;sumDeviationH = -50;
        cellP = cell(nBases,1);cellQ = cell(nBases,1);cellB = cell(nBases,1);
        cellM = cell(nBases,1);cellD = cell(nBases,1);cellBH = cell(nBases,1);W = cell(nUsers,1);
        
        xIteration = 0;
        scaContinue = 1;
                
        while scaContinue
            
            yIteration = 0;
            masterContinue = 1;
            
            if xIteration == 0
                for iBand = 1:nBands
                    for iUser = 1:nUsers
                        W{iUser,iBand} = ones(SimParams.nRxAntenna,nLayers) / sqrt(SimParams.nRxAntenna);
                    end
                end
            else
                for iBand = 1:nBands
                    for iUser = 1:nUsers
                        baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                        for iLayer = 1:nLayers
                            R = W{iUser,iBand}(:,iLayer) * W{iUser,iBand}(:,iLayer)' * SimParams.N;
                            for iBase = 1:nBases
                                for jUser = 1:usersPerCell(iBase,1)
                                    H = cH{iBase,iBand}(:,:,iUser);
                                    M = cellM{iBase,1}(:,:,jUser,iBand);
                                    R = R + H * (M * M') * H';
                                end
                            end
                            H = cH{baseNode,iBand}(:,:,iUser);
                            xUser = iUser == cellUserIndices{baseNode,1};
                            W{iUser,iBand}(:,iLayer) = R \ (H * cellM{baseNode,iBand}(:,iLayer,xUser,iBand));
                            W{iUser,iBand}(:,iLayer) = W{iUser,iBand}(:,iLayer) / norm(W{iUser,iBand}(:,iLayer),2);
                        end
                    end
                end
            end
            
            if xIteration == 0
                for iBase = 1:nBases
                    cellP{iBase,1} = ones(nLayers,usersPerCell(iBase,1),nBands) / usersPerCell(iBase,1);
                    cellQ{iBase,1} = ones(nLayers,usersPerCell(iBase,1),nBands) / usersPerCell(iBase,1);
                    cellB{iBase,1} = ones(nLayers,usersPerCell(iBase,1),nBands) * 10 + rand(nLayers,usersPerCell(iBase,1),nBands);
                end
            else
                for iBase = 1:nBases
                    cellB{iBase,1} = cellBH{iBase,1};
                    for iBand = 1:nBands
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:nLayers
                                cellP{iBase,1}(iLayer,iUser,iBand) = real(W{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * cellM{iBase,1}(:,iLayer,iUser,iBand));
                                cellQ{iBase,1}(iLayer,iUser,iBand) = imag(W{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * cellM{iBase,1}(:,iLayer,iUser,iBand));
                            end
                        end
                    end
                end
            end
            
            xIteration = xIteration + 1;
            currentIF = ones(nLayers,nUsers,nBases,nBands) * SimParams.sPower;
            
            if xIteration > mIterationsSCA
                scaContinue = 0;
            end
            
            while masterContinue
                
                yIteration = yIteration + 1;
                if yIteration > mIterationsSG
                    masterContinue = 0;
                end
                
                for iBase = 1:nBases
                    
                    kUsers = usersPerCell(iBase,1);
                    
                    cvx_begin
                    
                    dual variables dualD{nLayers,nUsers,nBands}
                    expressions p(nLayers,kUsers,nBands) q(nLayers,kUsers,nBands)
                    variable M(SimParams.nTxAntenna,nLayers,kUsers,nBands) complex
                    variables t(nLayers,kUsers,nBands) b(nLayers,kUsers,nBands) g(nLayers,kUsers,nBands)
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
                            H = cH{iBase,iBand}(:,:,cUser);
                            
                            for iLayer = 1:nLayers
                                
                                intVector = sqrt(SimParams.N) * W{cUser,iBand}(:,iLayer);
                                
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        intVector = [intVector ; sqrt(currentIF(iLayer,cUser,jBase,iBand))];
                                    end
                                end
                                
                                for jUser = 1:kUsers
                                    if jUser ~= iUser
                                        for jLayer = 1:nLayers
                                            intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                        end
                                    else
                                        for jLayer = 1:nLayers
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                            end
                                        end
                                    end
                                end
                                
                                dualD{iLayer,cUser,iBand} : norm(intVector,2) <= sqrt(b(iLayer,iUser,iBand));
                                log(1 + g(iLayer,iUser,iBand)) >= t(iLayer,iUser,iBand);
                                
                                p(iLayer,iUser,iBand) = real(W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                q(iLayer,iUser,iBand) = imag(W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                q(iLayer,iUser,iBand) == 0;
                                
                                (cellP{iBase,1}(iLayer,iUser,iBand)^2 + cellQ{iBase,1}(iLayer,iUser,iBand)^2) / (cellB{iBase,1}(iLayer,iUser,iBand)) + ...
                                    (2 / cellB{iBase,1}(iLayer,iUser,iBand)) * (cellP{iBase,1}(iLayer,iUser,iBand) * (p(iLayer,iUser,iBand) - cellP{iBase,1}(iLayer,iUser,iBand))) + ...
                                    (2 / cellB{iBase,1}(iLayer,iUser,iBand)) * (cellQ{iBase,1}(iLayer,iUser,iBand) * (q(iLayer,iUser,iBand) - cellQ{iBase,1}(iLayer,iUser,iBand))) - ...
                                    (cellP{iBase,1}(iLayer,iUser,iBand)^2 + cellQ{iBase,1}(iLayer,iUser,iBand)^2) / (2 * cellB{iBase,1}(iLayer,iUser,iBand)^2) * ...
                                    (b(iLayer,iUser,iBand) - cellB{iBase,1}(iLayer,iUser,iBand)) >= g(iLayer,iUser,iBand);
                                
                            end
                            
                        end
                        
                        norm(vec(M(:,:,:,iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                        
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                            for iLayer = 1:nLayers
                                intVector = [];
                                H = cH{iBase,iBand}(:,:,cUser);
                                for jUser = 1:kUsers
                                    for jLayer = 1:nLayers
                                        intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                    end
                                end
                                dualD{iLayer,cUser,iBand} : norm(intVector,2) <= sqrt(currentIF(iLayer,cUser,iBase,iBand));
                            end
                        end
                        
                    end
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        sum(vec(t(:,iUser,:))) <= QueuedPkts(cUser,1) * log(2);
                    end
                    
                    cvx_end
                    
                    if strfind(cvx_status,'Solved')
                        
                        cellM{iBase,1} = M;
                        cellBH{iBase,1} = b;
                        cellD{iBase,1} = dualD;
                        
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(:,iUser,:))) * log2(exp(1)),0);
                            SimParams.Debug.tempResource{2,1}{cUser,1} = [SimParams.Debug.tempResource{2,1}{cUser,1} sum(vec(t(:,iUser,:)))];
                            SimParams.Debug.tempResource{3,1}{cUser,1} = [SimParams.Debug.tempResource{3,1}{cUser,1} qDeviation];
                            for iBand = 1:nBands
                                SimParams.Debug.tempResource{4,1}{cUser,iBand} = [SimParams.Debug.tempResource{4,1}{cUser,iBand} sum(t(:,iUser,iBand))];
                            end
                        end
                        
                    end
                    
                end
                
                currentIFH = currentIF;
                currentIF = zeros(nLayers,nUsers,nBases,nBands);
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            for iLayer = 1:nLayers
                                cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                                baseNode = SimStructs.userStruct{cUser,1}.baseNode;
                                currentIF(iLayer,cUser,iBase,iBand) = currentIFH(iLayer,cUser,iBase,iBand) - alpha * (cellD{baseNode,1}{iLayer,cUser,iBand} - cellD{iBase,1}{iLayer,cUser,iBand});
                            end
                        end
                    end
                end
                
                currentIF = max(currentIF,0);
                if norm(vec(currentIF - currentIFH),2) <= 1e-3
                    masterContinue = 0;
                end
                
            end
            
            sumDeviation = sum(cell2mat(SimParams.Debug.tempResource{3,1}));
            if abs(sumDeviation(1,end) - sumDeviationH) < 1e-4
                scaContinue = 0;
            else
                sumDeviationH = sumDeviation(1,end);
            end
            
        end
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = cellM{iBase,1}(:,:,:,iBand);
            end
        end
        
    case 'DualMethod'
        
        alpha = 0.5;
        nLayers = SimParams.maxRank;
        mIterationsSCA = 50;mIterationsSG = 5;sumDeviationH = -50;
        cellP = cell(nBases,1);cellQ = cell(nBases,1);cellB = cell(nBases,1);
        cellM = cell(nBases,1);cellX = cell(nBases,1);cellBH = cell(nBases,1);W = cell(nUsers,1);
        
        xIteration = 0;
        scaContinue = 1;
                
        while scaContinue
            
            yIteration = 0;
            masterContinue = 1;
            if xIteration == 0
                for iBase = 1:nBases
                    cellP{iBase,1} = ones(nLayers,usersPerCell(iBase,1),nBands) / usersPerCell(iBase,1);
                    cellQ{iBase,1} = ones(nLayers,usersPerCell(iBase,1),nBands) / usersPerCell(iBase,1);
                    cellB{iBase,1} = ones(nLayers,usersPerCell(iBase,1),nBands) * 10 + rand(nLayers,usersPerCell(iBase,1),nBands);
                end
            else
                for iBase = 1:nBases
                    cellB{iBase,1} = cellBH{iBase,1};
                    for iBand = 1:nBands
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:nLayers
                                cellP{iBase,1}(iLayer,iUser,iBand) = real(cH{iBase,iBand}(:,:,cUser) * cellM{iBase,1}(:,iLayer,iUser,iBand));
                                cellQ{iBase,1}(iLayer,iUser,iBand) = imag(cH{iBase,iBand}(:,:,cUser) * cellM{iBase,1}(:,iLayer,iUser,iBand));
                            end
                        end
                    end
                end
            end
            
            if xIteration == 0
                for iBand = 1:nBands
                    for iUser = 1:nUsers
                        W{iUser,iBand} = ones(SimParams.nRxAntenna,nLayers) / sqrt(SimParams.nRxAntenna);
                    end
                end
            else
                for iBand = 1:nBands
                    for iUser = 1:nUsers
                        baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                        for iLayer = 1:nLayers
                            R = W{iUser,iBand}(:,iLayer) * W{iUser,iBand}(:,iLayer)' * SimParams.N;
                            for iBase = 1:nBases
                                for jUser = 1:usersPerCell(iBase,1)
                                    H = cH{iBase,iBand}(:,:,iUser);
                                    M = cellM{iBase,1}(:,:,jUser,iBand);
                                    R = R + H * (M * M') * H';
                                end
                            end
                            H = cH{baseNode,iBand}(:,:,iUser);
                            xUser = iUser == cellUserIndices{baseNode,1};
                            W{iUser,iBand}(:,iLayer) = R \ (H * cellM{baseNode,iBand}(:,iLayer,xUser,iBand));
                            W{iUser,iBand}(:,iLayer) = W{iUser,iBand}(:,iLayer) / norm(W{iUser,iBand}(:,iLayer),2);
                        end
                    end
                end
            end
            
            xIteration = xIteration + 1;
            currentDual = zeros(nLayers,nUsers,nBases,nBands);
            
            if xIteration > mIterationsSCA
                scaContinue = 0;
            end
            
            while masterContinue
                
                yIteration = yIteration + 1;
                if yIteration > mIterationsSG
                    masterContinue = 0;
                end
                
                for iBase = 1:nBases
                    
                    kUsers = usersPerCell(iBase,1);
                    
                    cvx_begin
                    
                    expressions p(nLayers,kUsers,nBands) q(nLayers,kUsers,nBands) tempFirst tempSecond
                    variable M(SimParams.nTxAntenna,nLayers,kUsers,nBands) complex
                    variables t(nLayers,kUsers,nBands) b(nLayers,kUsers,nBands) g(nLayers,kUsers,nBands) x(nLayers,nUsers,nBases,nBands)
                    variables userObjective(kUsers,1) epiObjective
                    
                    minimize(epiObjective)
                    
                    subject to
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        abs(QueuedPkts(cUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1);
                    end
                    
                    for iBand = 1:nBands
                        tempFirst = 0;
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for jBase = 1:nBases
                                if jBase ~= iBase
                                    tempFirst = tempFirst + sum(currentDual(:,cUser,jBase,iBand) .* x(:,cUser,jBase,iBand));
                                end
                            end                            
                        end 
                                      
                        tempSecond = 0;
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                            tempSecond = tempSecond + sum(currentDual(:,cUser,iBase,iBand) .* x(:,cUser,iBase,iBand));
                        end
                    end
                    
                    epiObjective >= norm(userObjective,1) + tempFirst + tempSecond;
                    
                    for iBand = 1:nBands
                        
                        for iUser = 1:kUsers
                            
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            H = cH{iBase,iBand}(:,:,cUser);
                            
                            for iLayer = 1:nLayers
                                
                                intVector = sqrt(SimParams.N) * W{cUser,iBand}(:,iLayer);
                                
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        intVector = [intVector ; x(iLayer,cUser,jBase,iBand)];
                                    end
                                end
                                
                                for jUser = 1:kUsers
                                    if jUser ~= iUser
                                        for jLayer = 1:nLayers
                                            intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                        end
                                    else
                                        for jLayer = 1:nLayers
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                            end
                                        end
                                    end
                                end
                                
                                norm(intVector,2) <= sqrt(b(iLayer,iUser,iBand));
                                log(1 + g(iLayer,iUser,iBand)) >= t(iLayer,iUser,iBand);
                                
                                p(iLayer,iUser,iBand) = real(W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                q(iLayer,iUser,iBand) = imag(W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                q(iLayer,iUser,iBand) == 0;
                                
                                (cellP{iBase,1}(iLayer,iUser,iBand)^2 + cellQ{iBase,1}(iLayer,iUser,iBand)^2) / (cellB{iBase,1}(iLayer,iUser,iBand)) + ...
                                    (2 / cellB{iBase,1}(iLayer,iUser,iBand)) * (cellP{iBase,1}(iLayer,iUser,iBand) * (p(iLayer,iUser,iBand) - cellP{iBase,1}(iLayer,iUser,iBand))) + ...
                                    (2 / cellB{iBase,1}(iLayer,iUser,iBand)) * (cellQ{iBase,1}(iLayer,iUser,iBand) * (q(iLayer,iUser,iBand) - cellQ{iBase,1}(iLayer,iUser,iBand))) - ...
                                    (cellP{iBase,1}(iLayer,iUser,iBand)^2 + cellQ{iBase,1}(iLayer,iUser,iBand)^2) / (2 * cellB{iBase,1}(iLayer,iUser,iBand)^2) * ...
                                    (b(iLayer,iUser,iBand) - cellB{iBase,1}(iLayer,iUser,iBand)) >= g(iLayer,iUser,iBand);
                                
                            end
                            
                        end
                        
                        norm(vec(M(:,:,:,iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                        
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                            for iLayer = 1:nLayers
                                intVector = [];
                                H = cH{iBase,iBand}(:,:,cUser);
                                for jUser = 1:kUsers
                                    for jLayer = 1:nLayers
                                        intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                    end
                                end
                                norm(intVector,2) <= x(iLayer,cUser,iBase,iBand);
                            end
                        end
                        
                    end
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        sum(vec(t(:,iUser,:))) <= QueuedPkts(cUser,1) * log(2);
                    end
                    
                    cvx_end
                    
                    if strfind(cvx_status,'Solved')
                        
                        cellM{iBase,1} = M;
                        cellBH{iBase,1} = b;
                        cellX{iBase,1} = x;
                        
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(:,iUser,:))) * log2(exp(1)),0);
                            SimParams.Debug.tempResource{2,1}{cUser,1} = [SimParams.Debug.tempResource{2,1}{cUser,1} sum(vec(t(:,iUser,:)))];
                            SimParams.Debug.tempResource{3,1}{cUser,1} = [SimParams.Debug.tempResource{3,1}{cUser,1} qDeviation];
                            for iBand = 1:nBands
                                SimParams.Debug.tempResource{4,1}{cUser,iBand} = [SimParams.Debug.tempResource{4,1}{cUser,iBand} sum(t(:,iUser,iBand))];
                            end
                        end
                        
                    end
                    
                end
                
                currentDualH = currentDual;
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:nLayers
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        currentDual(iLayer,cUser,jBase,iBand) = currentDualH(iLayer,cUser,jBase,iBand) - alpha * ...
                                            (cellX{iBase,1}(iLayer,cUser,jBase,iBand) - cellX{jBase,1}(iLayer,cUser,jBase,iBand));
                                    end
                                end
                            end
                        end
                    end
                end
                                            
                if norm(vec(currentDual - currentDualH),2) <= 1e-3
                    masterContinue = 0;
                end
                
            end
            
            sumDeviation = sum(cell2mat(SimParams.Debug.tempResource{3,1}));
            if abs(sumDeviation(1,end) - sumDeviationH) < 1e-4
                scaContinue = 0;
            else
                sumDeviationH = sumDeviation(1,end);
            end
            
        end
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = cellM{iBase,1}(:,:,:,iBand);
            end
        end
        
        
end

for iUser = 1:nUsers
    SimParams.Debug.tempResource{2,1}{iUser,1} = SimParams.Debug.tempResource{2,1}{iUser,1} * log2(exp(1));
    for iBand = 1:nBands
        SimParams.Debug.tempResource{4,1}{iUser,iBand} = SimParams.Debug.tempResource{4,1}{iUser,iBand} * log2(exp(1));
    end
end
