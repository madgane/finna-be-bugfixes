
function [SimParams,SimStructs] = getQWtdWSRM(SimParams,SimStructs)

nSB = SimParams.nBands;
cH = SimStructs.linkChan;

for iBase = 1:SimParams.nBases
    
    currentUsers = [];
    updatePrecoders = 'true';
    for iBand = 1:SimParams.nBands
        currentUsers = [currentUsers ; SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}];
    end
    
    currentUsers = unique(currentUsers);
    nUsers = length(currentUsers);
    Queues = zeros(nUsers,1);
    
    for iUser = 1:nUsers
        cUser = currentUsers(iUser,1);
        Queues(iUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
    end
    
    switch SimParams.weightedSumRateMethod
            
        case 'BandWMMSEApproach'
            
            updatePrecoders = 'false';
            currentDesign = SimParams.PrecodingMethod;
            currentApproach = SimParams.weightedSumRateMethod;
            
            SimParams.PrecodingMethod = 'Best_WMMSE_Method';
            SimParams.weightedSumRateMethod = 'PerformScheduling';
            
            for iBand = 1:SimParams.nBands
                
                for iUser = 1:nUsers
                    cUser = currentUsers(iUser,1);
                    SimStructs.userStruct{cUser,1}.weighingFactor = Queues(iUser,1);
                end
                
                SimParams.Debug.privateExchanges.takeOverBand = iBand;
                [SimParams, SimStructs] = getWeightedMMSEDesign(SimParams,SimStructs);
                
                [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,iBand);
                Queues = max(Queues - SimParams.Debug.privateExchanges.resAllocation(iBand,:)',0);
                
            end
            
            SimParams.PrecodingMethod = currentDesign;
            SimParams.weightedSumRateMethod = currentApproach;
            SimParams.privateExchanges = rmfield(SimParams.Debug.privateExchanges,'takeOverBand');
            
        case 'BalancingApproach'
            
            re_iterate = 1;cvx_hist = 100;
            qWeights = Queues / log2(exp(1));
            m_k_n_o = -rand(nUsers,nSB) * 10;
            
            while re_iterate
                
                cvx_begin
                
                expression m_k_n(nUsers,nSB)
                variable M_k_n(SimParams.nTxAntenna,nUsers,nSB) complex
                variables t_k_n(nUsers,nSB) b_k_n(nUsers,nSB) g_k_n(nUsers,nSB)
                variables obj_var(nUsers,1) obj_func_var
                
                minimize(obj_func_var)
                
                subject to
                
                for iUser = 1:nUsers
                    Queues(iUser,1) - sum(t_k_n(iUser,:)) * qWeights(iUser,1) <= obj_var(iUser,1);
                end
                
                obj_func_var >= norm(obj_var,2);
                
                for iUser = 1:nUsers
                    cUser = currentUsers(iUser,1);
                    for iSB = 1:nSB
                        if_vector = sqrt(SimParams.N);
                        cCH = cH{iBase,iSB}(:,:,cUser);
                        for jUser = 1:nUsers
                            if jUser ~= iUser
                                if_vector = [if_vector ; cCH * M_k_n(:,jUser,iSB)];
                            end
                        end
                        norm(if_vector,2) <= sqrt(b_k_n(iUser,iSB));
                        log(1 + g_k_n(iUser,iSB)) >= t_k_n(iUser,iSB) * qWeights(iUser,1);
                        
                        m_k_n(iUser,iSB) = g_k_n(iUser,iSB) - b_k_n(iUser,iSB);
                        4 * real(cCH * M_k_n(:,iUser,iSB)) + m_k_n_o(iUser,iSB)^2 + 2 * m_k_n_o(iUser,iSB) * (m_k_n(iUser,iSB) - m_k_n_o(iUser,iSB)) ...
                            >= (g_k_n(iUser,iSB) + b_k_n(iUser,iSB))^2;
                        imag(cCH * M_k_n(:,iUser,iSB)) == 0;
                    end
                    norm(t_k_n(iUser,:),1) <= 1;
                end
                
                for iSB = 1:nSB
                    norm(vec(M_k_n(:,:,iSB)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iSB));
                end
                
                g_k_n >= 0;
                t_k_n >= 0;
                t_k_n <= 1;
                
                cvx_end
                
                if strcmp(cvx_status,'Solved')
                    m_k_n_o = g_k_n - b_k_n;
                    if abs(cvx_optval - cvx_hist) <= 1e-5
                        re_iterate = 0;
                    else
                        cvx_hist = cvx_optval;
                    end
                else
                    m_k_n_o = m_k_n_o * 2;
                end
                
            end
            
        otherwise
            
            display('Undefined Optimization Approach !');
            
    end
    
    if strcmp(updatePrecoders,'true')
        for iSB = 1:nSB
            SimStructs.baseStruct{iBase,1}.P{iSB,1} = M_k_n(:,:,iSB);
        end
    end            
    
end

