
function [SimParams,SimStructs] = getQWtdWSRM(SimParams,SimStructs)

nSB = SimParams.nBands;
cH = SimStructs.linkChan;

for iBase = 1:SimParams.nBases
    
    currentUsers = [];
    for iBand = 1:SimParams.nBands
        currentUsers = [currentUsers ; SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}];
    end
    
    iter_count = 0;
    max_iterations = 100;
    currentUsers = unique(currentUsers);
    nUsers = length(currentUsers);
    Queues = zeros(nUsers,1);
    
    for iUser = 1:nUsers
        cUser = currentUsers(iUser,1);
        Queues(iUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
    end
    
    switch SimParams.weightedSumRateMethod
        
        case 'GMApproach'
            
            th_tol = 1e-3;
            t_k_hist = 1e5;
            re_iterate = 1;
            phi = rand(nUsers,nSB);
            
            while re_iterate
                
                cvx_begin
                
                variables min_var z_var
                variables t_k(nUsers,1) mu_k_n(nUsers,nSB) beta_k_n(nUsers,nSB)
                variable M_k_n(SimParams.nTxAntenna,nUsers,nSB) complex;
                expressions x_1 x_2 beta_IF
                
                dual y1
                
                minimize(min_var - log(z_var) * nSB * nUsers * log2(exp(1)))
                
                subject to
                
                norm(t_k,2) <= min_var;
                y1:geo_mean(vec(mu_k_n)) >= z_var;
                
                for iSB = 1:nSB
                    norm(vec(M_k_n(:,:,iSB)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iSB));
                end
                
                for iUser = 1:nUsers
                    geo_mean(mu_k_n(iUser,:)) >= 2^((Queues(iUser,1) - t_k(iUser,1)) / nSB);
                end
                
                for iUser = 1:nUsers
                    icUser = currentUsers(iUser,1);
                    for iSB = 1:nSB
                        betaIF = sqrt(SimParams.N);
                        for jUser = 1:nUsers
                            if jUser ~= iUser
                                betaIF = [betaIF ; cH{iBase,iSB}(:,:,icUser) * M_k_n(:,jUser,iSB)];
                            end
                        end
                        
                        norm(betaIF,2) <= beta_k_n(iUser,iSB);
                        
                        x_1 = (phi(iUser,iSB) / 2) * beta_k_n(iUser,iSB)^2;
                        x_2 = (1 / (2 * phi(iUser,iSB))) * (mu_k_n(iUser,iSB) - 1);
                        
                        x_1 + x_2  <= real(cH{iBase,iSB}(:,:,icUser) * M_k_n(:,iUser,iSB));
                        imag(cH{iBase,iSB}(:,:,icUser) * M_k_n(:,iUser,iSB)) == 0;
                        
                        mu_k_n(iUser,iSB) >= 1;
                        
                    end
                end
                
                cvx_end
                
                [t_k mu_k_n Queues  sum(log2(mu_k_n),2)]
                
                if strfind(cvx_status,'Solved')
                    for iUser = 1:nUsers
                        for iSB = 1:nSB
                            phi(iUser,iSB) = real(sqrt(mu_k_n(iUser,iSB) - 1) / beta_k_n(iUser,iSB));
                        end
                    end
                    
                    iter_count = iter_count + 1;
                    if abs(norm(t_k,2) - t_k_hist) < th_tol
                        re_iterate = 0;
                    end
                    t_k_hist = norm(t_k,2);
                else
                    phi = phi./2;
                end
                
                if iter_count > max_iterations
                    re_iterate = 0;
                end
                
            end
            
            for iSB = 1:nSB
                SimStructs.baseStruct{iBase,1}.P{iSB,1} = M_k_n(:,:,iSB);
            end
            
        case 'DCApproach'
            
            th_tol = 1e-3;
            t_k_hist = 1e5;
            re_iterate = 1;
            
            i_ur_k_n = zeros(nUsers,nSB);
            i_ui_k_n = zeros(nUsers,nSB);
            i_beta_k_n = zeros(nUsers,nSB);
            
            for iUser = 1:nUsers
                icUser = currentUsers(iUser,1);
                for iSB = 1:nSB
                    i_ur_k_n(iUser,iSB) = norm(cH{iBase,iSB}(:,:,icUser))^2;
                    i_beta_k_n(iUser,iSB) = SimParams.N;
                    for jUser = 1:nUsers
                        jcUser = currentUsers(jUser,1);
                        if jUser ~= iUser
                            i_beta_k_n(iUser,iSB) = i_beta_k_n(iUser,iSB) + norm(cH{iBase,iSB}(:,:,jcUser))^2 / nUsers;
                        end
                    end
                end
            end
            
            phi = i_ur_k_n ./ i_beta_k_n;
            
            cvx_quiet('false');
            
            while re_iterate
                
                cvx_begin
                
                variables min_var z_var
                variables t_k(nUsers,1) mu_k_n(nUsers,nSB) beta_k_n(nUsers,nSB)
                variable M_k_n(SimParams.nTxAntenna,nUsers,nSB) complex
                expressions beta_IF ur_k_n(nUsers,nSB) ui_k_n(nUsers,nSB)

                minimize(min_var - log(z_var) * nSB * nUsers * log2(exp(1)))
                
                subject to
                
                norm(t_k,2) <= min_var;
                z_var <= geo_mean(vec(mu_k_n));
                
                for iSB = 1:nSB
                    norm(vec(M_k_n(:,:,iSB)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iSB));
                end
                
                for iUser = 1:nUsers
                    geo_mean(mu_k_n(iUser,:)) >= 2^((Queues(iUser,1) - t_k(iUser,1)) / nSB);
                end
                
                for iUser = 1:nUsers
                    icUser = currentUsers(iUser,1);
                    for iSB = 1:nSB
                        betaIF = sqrt(SimParams.N);
                        for jUser = 1:nUsers
                            if jUser ~= iUser
                                betaIF = [betaIF ; cH{iBase,iSB}(:,:,icUser) * M_k_n(:,jUser,iSB)];
                            end
                        end
                        
                        norm(betaIF,2) <= sqrt(beta_k_n(iUser,iSB));
                        
                        ur_k_n(iUser,iSB) = real(cH{iBase,iSB}(:,:,icUser) * M_k_n(:,iUser,iSB));
                        ui_k_n(iUser,iSB) = imag(cH{iBase,iSB}(:,:,icUser) * M_k_n(:,iUser,iSB));
                        
                        i_ur_k_n(iUser,iSB)^2 + i_ui_k_n(iUser,iSB)^2 + i_beta_k_n(iUser,iSB) ...
                            + 2 * (i_ur_k_n(iUser,iSB) * (ur_k_n(iUser,iSB) - i_ur_k_n(iUser,iSB)) ...
                            + i_ui_k_n(iUser,iSB) * (ui_k_n(iUser,iSB) - i_ui_k_n(iUser,iSB))) ...
                            + (beta_k_n(iUser,iSB) - i_beta_k_n(iUser,iSB)) ...
                            >= mu_k_n(iUser,iSB)^2 / (2 * phi(iUser,iSB)) + (phi(iUser,iSB) / 2) * beta_k_n(iUser,iSB)^2;
                        
                        imag(cH{iBase,iSB}(:,:,icUser) * M_k_n(:,iUser,iSB)) == 0;
                        mu_k_n(iUser,iSB) >= 1;
                        
                    end
                end
                
                cvx_end
                
                [t_k mu_k_n Queues  sum(log2(mu_k_n),2)]
                
                if strfind(cvx_status,'Solved')
                    
                    i_ur_k_n = ur_k_n;
                    i_ui_k_n = ui_k_n;
                    i_beta_k_n = beta_k_n;
                    phi = mu_k_n ./ beta_k_n;
                    
                    iter_count = iter_count + 1;
                    if abs(norm(t_k,2) - t_k_hist) < th_tol
                        re_iterate = 0;
                    end
                    t_k_hist = norm(t_k,2);
                    
                else
                    phi = phi./2;
                    i_beta_k_n = i_beta_k_n * 10;
                    i_ur_k_n = i_ur_k_n / 2;
                end
                
                if iter_count > max_iterations
                    re_iterate = 0;
                end
                
            end
            
            for iSB = 1:nSB
                SimStructs.baseStruct{iBase,1}.P{iSB,1} = M_k_n(:,:,iSB);
            end
            
        case 'BandWMMSEApproach'
            
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
                [SimParams SimStructs] = getWeightedMMSEDesign(SimParams,SimStructs);
                
                [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,iBand);
                Queues = max(Queues - SimParams.Debug.privateExchanges.resAllocation(iBand,:)',0);
                
            end
            
            SimParams.PrecodingMethod = currentDesign;
            SimParams.weightedSumRateMethod = currentApproach;
            SimParams.privateExchanges = rmfield(SimParams.Debug.privateExchanges,'takeOverBand');
            
        case 'JointSchedWMMSEApproach'
            
            xChannel = zeros(SimParams.nTxAntenna,length(currentUsers));
            
            for iBand = 1:SimParams.nBands
                
                for iUser = 1:nUsers
                    cUser = currentUsers(iUser,1);
                    xChannel(:,iUser) = Queues(iUser,1) * cH{iBase,iBand}(:,:,cUser).';
                    SimStructs.userStruct{cUser,1}.weighingFactor = Queues(iUser,1);
                end
                
                SimParams.Debug.privateExchanges.takeOverBand = iBand;
                [SimParams SimStructs] = getStrWeightedMMSEDesign(SimParams,SimStructs);
                
                [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,iBand);
                Queues = max(Queues - SimParams.Debug.privateExchanges.resAllocation(iBand,:)',0);
                
            end
            
            SimParams.PrecodingMethod = currentDesign;
            SimParams.weightedSumRateMethod = currentApproach;
            SimParams.privateExchanges = rmfield(SimParams.Debug.privateExchanges,'takeOverBand');
            
        case 'AltWMMSEApproach'
            
            th_tol = 1e-3;
            t_k_hist = 1e5;
            re_iterate = 1;
            priorityWts = ones(nUsers,nSB) / nSB;
            
            while re_iterate
                
                for iBand = 1:SimParams.nBands
                    
                    for iUser = 1:nUsers
                        cUser = currentUsers(iUser,1);
                        SimStructs.userStruct{cUser,1}.weighingFactor = priorityWts(iUser,iBand);
                    end
                    
                    SimParams.Debug.privateExchanges.takeOverBand = iBand;
                    [SimParams SimStructs] = getWeightedMMSEDesign(SimParams,SimStructs);
                    [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,iBand);
                    
                end
                
                rateAchieved = SimParams.Debug.privateExchanges.resAllocation;
                
                cvx_begin
                
                variables wts_k_n(nUsers,nSB) t_k(nUsers,1)
                
                minimize(norm((Queues - t_k),2))
                
                subject to
                
                for iUser = 1:nUsers
                    wts_k_n(iUser,:) * rateAchieved(:,iUser) >= t_k(iUser,1);
                    norm(wts_k_n(iUser,:),1) <= 1;
                end
                
                wts_k_n >= 0;
                
                cvx_end
                
                priorityWts = wts_k_n;
                
                if abs(norm(t_k,2) - t_k_hist) / norm(t_k,2) < th_tol
                    re_iterate = 0;
                else
                    t_k_hist = norm(t_k,2);
                end
                
            end
            
            SimParams.privateExchanges = rmfield(SimParams.Debug.privateExchanges,'takeOverBand');
            
        case 'SINRBalancing'
            
            th_tol = 1e-3;
            t_k_hist = 1e5;
            re_iterate = 1;
            phi = rand(nUsers,nSB);
            
            while re_iterate
                
                cvx_begin

                variable min_var
                variables t_k(nUsers,1) mu_k_n(nUsers,nSB) beta_k_n(nUsers,nSB)
                variable M_k_n(SimParams.nTxAntenna,nUsers,nSB) complex
                expressions beta_IF ur_k_n(nUsers,nSB) ui_k_n(nUsers,nSB)
                
                minimize(min_var)
                
                subject to
                
                norm(Queues - t_k .* (2.^(Queues ./ nSB)),2) <= min_var;
                
                for iSB = 1:nSB
                    norm(vec(M_k_n(:,:,iSB)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iSB));
                end
                
                for iUser = 1:nUsers
                    geo_mean(mu_k_n(iUser,:)) >= t_k(iUser,1) * 2^(Queues(iUser,1) / nSB);
                end
                
                for iUser = 1:nUsers
                    icUser = currentUsers(iUser,1);
                    for iSB = 1:nSB
                        betaIF = sqrt(SimParams.N);
                        for jUser = 1:nUsers
                            if jUser ~= iUser
                                betaIF = [betaIF ; cH{iBase,iSB}(:,:,icUser) * M_k_n(:,jUser,iSB)];
                            end
                        end
                        
                        norm(betaIF,2) <= beta_k_n(iUser,iSB);
                        x_1 = (phi(iUser,iSB) / 2) * beta_k_n(iUser,iSB)^2;
                        x_2 = (1 / (2 * phi(iUser,iSB))) * (mu_k_n(iUser,iSB) - 1);
                        
                        x_1 + x_2  <= real(cH{iBase,iSB}(:,:,icUser) * M_k_n(:,iUser,iSB));
                        imag(cH{iBase,iSB}(:,:,icUser) * M_k_n(:,iUser,iSB)) == 0;
                        
                        mu_k_n(iUser,iSB) >= 1;
                        
                    end
                end
                
                t_k >= 0;
                t_k <= 1;
                
                cvx_end
                
                [t_k mu_k_n Queues  sum(log2(mu_k_n),2)]
                
                if strfind(cvx_status,'Solved')
                    for iUser = 1:nUsers
                        for iSB = 1:nSB
                            phi(iUser,iSB) = sqrt(mu_k_n(iUser,iSB) - 1) / beta_k_n(iUser,iSB);
                        end
                    end
                    
                    iter_count = iter_count + 1;
                    if abs(min_var - t_k_hist) < th_tol
                        re_iterate = 0;
                    end
                    t_k_hist = min_var;
                else
                    phi = phi./2;
                end
                
                if iter_count > max_iterations
                    re_iterate = 0;
                end
                
            end
            
            for iSB = 1:nSB
                SimStructs.baseStruct{iBase,1}.P{iSB,1} = M_k_n(:,:,iSB);
            end

            
        otherwise
            
            display('Undefined Optimization Approach !');
            
    end
    
end

