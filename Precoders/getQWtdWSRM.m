
function [SimParams,SimStructs] = getQWtdWSRM(SimParams,SimStructs)

nSB = SimParams.nBands;
cH = SimStructs.linkChan;

for iBase = 1:SimParams.nBases
    
    currentUsers = [];
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
    
    th_tol = 1e-3;
    t_k_hist = 1e5;
    re_iterate = 1;
    phi = rand(nUsers,nSB);
    
    while re_iterate
        
        cvx_begin
        
        variable min_var
        variables t_k(nUsers,1) mu_k_n(nUsers,nSB) beta_k_n(nUsers,nSB)
        variable M_k_n(SimParams.nTxAntenna,nUsers,nSB) complex
        expressions x_1 x_2 beta_IF
        
        minimize(min_var)
        
        subject to
        
        norm(t_k,2) <= min_var;
        
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
                    phi(iUser,iSB) = sqrt(mu_k_n(iUser,iSB) - 1) / beta_k_n(iUser,iSB);
                end
            end
            
            if abs(norm(t_k,2) - t_k_hist) < th_tol
                re_iterate = 0;
            end            
            t_k_hist = norm(t_k,2);
        else            
            phi = phi./2;            
        end
        
    end    
    
    for iSB = 1:nSB
        SimStructs.baseStruct{iBase,1}.P{iSB,1} = M_k_n(:,:,iSB);
    end
    
end

end
