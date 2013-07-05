
function [SimParams SimStructs] = getWSRM_SCP(SimParams,SimStructs)

for iBand = 1:SimParams.nBands
    
    kUsers = 0;
    txUserSet = cell(SimParams.nBases,1);
    txStreamSet = cell(SimParams.nBases,1);
    
    for iBase = 1:SimParams.nBases
        txUserSet{iBase,1} = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1};
        txStreamSet{iBase,1} = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1};
        kUsers = kUsers + length(txUserSet{iBase,1});
    end
    
    H = SimStructs.linkChan;
    kUserStreams = kUsers * SimParams.maxRank;
    
    alpha = ones(kUserStreams,1) * 1.5;
    tOld = ones(kUserStreams,1) + rand(kUserStreams,1);
    phi = sqrt(tOld - 1) ./ SimParams.N;
    
    scalingfactor = 1.5;
    while(min(alpha) <= 1)
        weightscaled = weightscaled * scalingfactor;
    end
   
    z_old = 0;
    max_iter = 100;
    tolerance = 1e-3;
    const_iter = 1;
    
    W = cell(kUserStreams,1);
    for iUser = 1:kUserStreams
        W{iUser,1} = ones(SimParams.nRxAntenna,1);
    end
    
    cvx_solver sdpt3
    cvx_quiet(true);
    
    while const_iter
        
        cvx_begin
        
        variable z
        variable X(SimParams.nTxAntenna,kUserStreams) complex
        variable t(kUserStreams,1)
        variable x(kUserStreams,1)
        variable beta_var(kUserStreams,1)
        
        maximize(z)
        
        subject to
        
        z <= geo_mean(t);
        
        for iBase = 1:SimParams.nBases
            B = [];
            for iUser = 1:length(txUserSet{iBase,1})
                
                iiUser = txUserSet{iBase,1}(iUser,1);
                
                for iStream = 1:SimParams.maxRank
                    
                    G = SimParams.N;
                    cIndex = (iBase - 1) * length(txUserSet{iBase,1}) * SimParams.maxRank + (iUser - 1) * SimParams.maxRank + iStream;
                    
                    S = real(W{cIndex,1}' * H{iBase,iBand}(:,:,iiUser) * X(:,cIndex));
                    vector = [0.5 * (S - (0.5 * x(cIndex,1) / phi(cIndex,1)) - 1) ; sqrt(0.5 * phi(cIndex,1)) * beta_var(cIndex,1)];
                    
                    {vector,0.5 * (S - (0.5 * x(cIndex,1) / phi(cIndex,1)) + 1)} <In> complex_lorentz(length(vector));
                    tOld(cIndex,1)^(1/alpha(cIndex,1)) + (1/alpha(cIndex,1)) * tOld(cIndex,1)^(1/alpha(cIndex,1) - 1) * (t(cIndex,1) - tOld(cIndex,1)) <= x(cIndex,1) + 1;
                    
                    for jBase = 1:SimParams.nBases
                        if jBase ~= iBase
                            for jUser = 1:length(txUserSet{jBase,1})
                                for jStream = 1:SimParams.maxRank
                                    kIndex = (jBase - 1) * length(txUserSet{jBase,1}) * SimParams.maxRank + ((jUser - 1) * SimParams.maxRank + jStream);
                                    S = W{cIndex,1}' * H{jBase,iBand}(:,:,iiUser) * X(:,kIndex);
                                    G = [G ; S];
                                end
                            end
                        end
                    end
                    
                    for jUser = 1:length(txUserSet{iBase,1})
                        if jUser ~= iUser
                            for jStream = 1:SimParams.maxRank
                                kIndex = (iBase - 1) * length(txUserSet{iBase,1}) * SimParams.maxRank + ((jUser - 1) * SimParams.maxRank + jStream);
                                S = W{cIndex,1}' * H{iBase,iBand}(:,:,iiUser) * X(:,kIndex);
                                G = [G ; S];
                            end
                        end
                    end
                    
                    {G,beta_var(cIndex,1)} <In> complex_lorentz(length(G));
                    B = [B ; X(:,cIndex,1)];
                    
                    imag(W{cIndex,1}' * H{iBase,iBand}(:,:,iiUser) * X(:,cIndex)) == 0;
                    
                end
            end
            
            norm(B) <= sqrt(SimParams.sPower);
            
        end
        
        t >= 1;
        z >= 0;
        x >= 0;
        
        cvx_end
        
        if strcmp(cvx_status,'Solved')
            tOld = t;
            phi = sqrt(x)./beta_var;

            if abs(z_old - z) < tolerance
                const_iter = 0;
            else
                z_old = z;
                const_iter = const_iter + 1;
            end
            
            if const_iter > max_iter
                const_iter = 0;
            end
            
            if strcmp(SimParams.plotMode,'network_rate_convergence')
                if const_iter >= SimParams.iDrop
                    const_iter = 0;
                end
            end

        else
            phi = phi / 2; 
            const_iter = 1;
        end
                      
    end
    
    iIndex = 0;
    for iBase = 1:SimParams.nBases
        P = [];
        for iUser = 1:length(txUserSet{iBase,1})
            Ww = [];
            for iStream = 1:SimParams.maxRank
                iIndex = iIndex + 1;
                P = [P X(:,iIndex)];
                Ww = [Ww W{iIndex,1}];
            end
            SimStructs.userStruct{txUserSet{iBase,1}(iUser,1)}.W = Ww;         
        end
        SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
    end
end

