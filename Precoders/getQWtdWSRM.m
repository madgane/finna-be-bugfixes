
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
                
                minimize(min_var)
                
                subject to
                
                norm(t_k,2) <= min_var;
                
                for iSB = 1:nSB
                    norm(vec(M_k_n(:,:,iSB)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iSB));
                end
                
                for iUser = 1:nUsers
                    geo_mean(mu_k_n(iUser,:)) >= (2^((Queues(iUser,1) - t_k(iUser,1))))^(1 / nSB);
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
                    geo_mean(mu_k_n(iUser,:)) >= (2^((Queues(iUser,1) - t_k(iUser,1))))^(1 / nSB);
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
            
        case 'QwtWSRMApproach'
            
            U_i_o = zeros(nUsers,nSB);
            U_r_o = rand(nUsers,nSB) * 2 + 1;
            beta_k_n_o = rand(nUsers,nSB) * 2 + 2;
            t_k_n_o = U_r_o ./ beta_k_n_o;
            
            re_iterate = 1;
            threshold = 1e-3;
            prod_history = 1e10;
            qWeights = Queues;
            
            while re_iterate
                
                cvx_begin
                
                expressions prod_t U_r(nUsers,nSB) U_i(nUsers,nSB)
                variables beta_k_n(nUsers,nSB) t_k_n(nUsers,nSB) max_var(nSB,1) max_obj
                variable M_k_n(SimParams.nTxAntenna,SimParams.nRxAntenna,nUsers,nSB) complex
                
                maximize(max_obj)
                
                subject to
                
                for iSB = 1:nSB
                    prod_t = t_k_n(:,iSB);
                    {prod_t,max_var(iSB,1)} == geo_mean_cone(nUsers,1,qWeights,'POS');
                end
                
                {max_var,max_obj} == geo_mean_cone(nSB,1,'POS');
                
                for iUser = 1:nUsers
                    cUser = currentUsers(iUser,1);
                    
                    for iSB = 1:nSB
                        if_vector = sqrt(SimParams.N);
                        currentH = cH{iBase,iSB}(:,:,cUser);
                        
                        for jUser = 1:nUsers
                            if jUser ~= iUser
                                if_vector = [if_vector ; currentH * M_k_n(:,:,jUser,iSB)];
                            end
                        end
                        
                        norm(if_vector,2) <= sqrt(beta_k_n(iUser,iSB));
                        U_r(iUser,iSB) = real(currentH * M_k_n(:,:,iUser,iSB));U_i(iUser,iSB) = imag(currentH * M_k_n(:,:,iUser,iSB));
                        
                        U_r_o(iUser,iSB)^2 + U_i_o(iUser,iSB)^2 + 2 * ...
                            (U_r_o(iUser,iSB) * (U_r(iUser,iSB) - U_r_o(iUser,iSB)) + U_i_o(iUser,iSB) * (U_i(iUser,iSB) - U_i_o(iUser,iSB))) ...
                            + beta_k_n(iUser,iSB) >= (t_k_n_o(iUser,iSB) * beta_k_n_o(iUser,iSB) / 2) * ...
                            ((t_k_n(iUser,iSB) / t_k_n_o(iUser,iSB))^2 + (beta_k_n(iUser,iSB) / beta_k_n_o(iUser,iSB))^2);
                        
                        U_i == 0;
                        
                    end
                    
                    prod(t_k_n_o(iUser,:).^nSB) * sum((1./nSB) .* pow_pos(t_k_n(iUser,:) ./ t_k_n_o(iUser,:),nSB)) <= 2^Queues(iUser,1);
                    
                end
                
                for iSB = 1:nSB
                    x_vec = [];
                    for iUser = 1:nUsers
                        x_vec = [x_vec ; vec(M_k_n(:,:,iUser,iSB))];
                    end
                    {x_vec,sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iSB))} == complex_lorentz(length(x_vec));
                end
                
                t_k_n >= 1;
                
                cvx_end
                
                [Queues sum(log2(t_k_n),2)]
                
                if strfind(cvx_status,'Solved')
                    t_k_n_o = abs(t_k_n);
                    beta_k_n_o = abs(beta_k_n);
                    U_r_o = U_r;U_i_o = U_i;
                    
                    if abs(prod_history - cvx_optval) <= threshold
                        re_iterate = 0;
                    else
                        prod_history = cvx_optval;
                    end
                    
                else
                    beta_k_n_o = beta_k_n_o * 2;
                    U_r_o = U_r_o / 2;
                    t_k_n_o = max(U_r_o ./ beta_k_n_o,1);
                end
                
            end
            
            for iSB = 1:nSB
                for iUser = 1:nUsers
                    SimStructs.baseStruct{iBase,1}.P{iSB,1} = [SimStructs.baseStruct{iBase,1}.P{iSB,1} M_k_n(:,:,iUser,iSB)];
                end
            end
            
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
%                 maximize(obj_func_var);
                
                subject to
                
                for iUser = 1:nUsers
                    Queues(iUser,1) - sum(t_k_n(iUser,:)) * qWeights(iUser,1) <= obj_var(iUser,1);
%                     qWeights(iUser,1) * sum(t_k_n(iUser,:)) >= obj_var(iUser,1);
                end
                
                obj_func_var >= norm(obj_var,2);
%                 obj_func_var <= sum(obj_var);
                
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
            
            for iSB = 1:nSB
                SimStructs.baseStruct{iBase,1}.P{iSB,1} = M_k_n(:,:,iSB);
            end
            
        otherwise
            
            display('Undefined Optimization Approach !');
            
    end
    
end

