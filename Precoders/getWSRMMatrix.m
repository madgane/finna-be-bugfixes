
function [SimParams SimStructs] = getWSRMMatrix(SimParams,SimStructs)

for iBand = 1:SimParams.nBands
    
    switch SimParams.weightedSumRateMethod
        
        case 'GMApproach'
            
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
            cvx_quiet(false);
            
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
                    
                    norm(B) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    
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
            
            
        case 'TaylorApproach'
            
            totUsers = zeros(SimParams.nBases,1);
            currentUsers = cell(SimParams.nBases,1);
            currentWeights = cell(SimParams.nBases,1);
            
            for iBase = 1:SimParams.nBases
                currentUsers{iBase,1} = unique(SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1});
                totUsers(iBase,1) = length(currentUsers{iBase,1});
                currentWeights{iBase,1} = zeros(length(currentUsers{iBase,1}),1);
                for iUser = 1:totUsers(iBase,1)
                    cUser = currentUsers{iBase,1}(iUser,1);
                    currentWeights{iBase,1}(iUser,1) = SimStructs.userStruct{cUser,1}.weighingFactor;
                end
            end
            
            re_iterate = 1;
            threshold = 1e-1;
            obj_var_hist = 10;
            cH = SimStructs.linkChan;
            
            t = cell(SimParams.nBases,1);
            M = cell(SimParams.nBases,1);
            beta = cell(SimParams.nBases,1);
            
            obj_var = sdpvar(1);
            for iBase = 1:SimParams.nBases
                t{iBase,1} = sdpvar(totUsers(iBase,1),1);
                beta{iBase,1} = sdpvar(totUsers(iBase,1),1);
                M{iBase,1} = sdpvar(SimParams.nTxAntenna,SimParams.nRxAntenna,totUsers(iBase,1),'full','complex');
            end
            
            t_o = cell(SimParams.nBases,1);
            beta_o = cell(SimParams.nBases,1);
            
            for iBase = 1:SimParams.nBases
                t_o{iBase,1} = 2 * rand(totUsers(iBase,1),1) + 1;
                beta_o{iBase,1} = rand(totUsers(iBase,1),1) + 10;
            end
            
            problem_options = sdpsettings('solver','sdpt3','verbose',0);
            
            while re_iterate
                
                current_objective = [];%#ok<*AGROW>
                problem_constraints = [];%#ok<*AGROW>
                
                for iBase = 1:SimParams.nBases
                    
                    for iUser = 1:totUsers(iBase,1)
                        
                        cUser = currentUsers{iBase,1}(iUser,1);
                        H = cH{iBase,iBand}(:,:,cUser);
                        
                        c_tko = t_o{iBase,1}(iUser,1);
                        c_bko = beta_o{iBase,1}(iUser,1);
                        c_iwt = 1 / currentWeights{iBase,1}(iUser,1);
                        
                        current_constraint = [real(H * M{iBase,1}(:,:,iUser)) >= ...
                            (0.5 / sqrt(c_tko^(c_iwt) - 1)) * c_iwt * (c_tko^(c_iwt - 1)) * c_bko * (t{iBase,1}(iUser,1) - c_tko) ...
                            + sqrt(c_tko^(c_iwt) - 1) * (beta{iBase,1}(iUser,1) - c_bko)];
                        
                        problem_constraints = problem_constraints + current_constraint;
                        problem_constraints = problem_constraints + [imag(H * M{iBase,1}(:,:,iUser)) == 0];
                        
                        sumIF = [SimParams.N];
                        for jBase = 1:SimParams.nBases
                            if jBase == iBase
                                for jUser = 1:totUsers(jBase,1)
                                    if jUser ~= iUser
                                        sumIF = [sumIF , H * M{jBase,1}(:,:,jUser)];
                                    end
                                end
                            else
                                Hif = cH{jBase,iBand}(:,:,cUser);
                                for jUser = 1:totUsers(jBase,1)
                                    sumIF = [sumIF , Hif * M{jBase,1}(:,:,jUser)];
                                end
                            end
                        end
                        
                        current_constraint = [cone(sumIF,beta{iBase,1}(iUser,1))];
                        problem_constraints = problem_constraints + current_constraint;
                        
                    end
                    
                    totPower = [];
                    for iUser = 1:totUsers(iBase,1)
                        totPower = [totPower ; vec(M{iBase,1}(:,:,iUser))];
                    end
                    
                    current_constraint = [cone(totPower,sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand)))];
                    problem_constraints = problem_constraints + current_constraint;
                    
                    current_objective = [current_objective ; t{iBase,1}];
                    problem_constraints = problem_constraints + [t{iBase,1} >= 1];
                    problem_constraints = problem_constraints + [beta{iBase,1} >= 1];
                    
                end
                
                problem_constraints = problem_constraints + [geomean(current_objective) >= obj_var];
                
                problem_objective = -obj_var;
                problem_solution = solvesdp(problem_constraints,problem_objective,problem_options);
                
                yalmiperror(problem_solution.problem)
                [log2(double(t{1,1})) double(beta{1,1})]
                
                switch problem_solution.problem
                    
                    case {0,3}
                        
                        for iBase = 1:SimParams.nBases
                            t_o{iBase,1} = double(t{iBase,1});
                            beta_o{iBase,1} = double(beta{iBase,1});
                        end
                        
                        if (abs(log2(double(obj_var)) - obj_var_hist)) <= threshold
                            re_iterate = 0;
                            for iBase = 1:SimParams.nBases
                                P = [];
                                for iUser = 1:totUsers(iBase,1)
                                    P = [P double(M{iBase,1}(:,:,iUser))];
                                end
                                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
                            end
                        else
                            obj_var_hist = log2(double(obj_var));
                        end
                        
                    otherwise
                        
                        for iBase = 1:SimParams.nBases
                            beta_o{iBase,1} = beta_o{iBase,1} * 2;
                        end
                end
                
            end
            
        case 'LousyApproach'
            
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
            cvx_quiet(false);
            
            while const_iter
                
                cvx_begin
                
                variable z
                variable X(SimParams.nTxAntenna,kUserStreams) complex
                variable t(kUserStreams,1)
                variable x(kUserStreams,1)
                variable beta_var(kUserStreams,1)
                
                maximize(z)
                
                subject to
                
                {t,z} == geo_mean_cone(kUserStreams,1,alpha,'POS');                
                
                for iBase = 1:SimParams.nBases
                    B = [];
                    for iUser = 1:length(txUserSet{iBase,1})
                        
                        iiUser = txUserSet{iBase,1}(iUser,1);
                        
                        for iStream = 1:SimParams.maxRank
                            

                            
                            real(W{cIndex,1}' * H{iBase,iBand}(:,:,iiUser) * X(:,cIndex))^2
                            
                            
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
                    
                    norm(B) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    
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
            
            
        otherwise
            
            display('Unknown WSRM Scheme !');
            
    end
    
end

