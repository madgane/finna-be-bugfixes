function [SumCapacity] = performMockReception(SimParams,SimStructs,V,iBand)

SumCapacity = 0;
H = SimStructs.linkChan;

for iUser = 1:SimParams.nUsers
    
    if ~isempty(V{iUser,1})
        
        cUser = SimStructs.userStruct{iUser,1};    
        N = SimParams.N * eye(SimParams.nRxAntenna);
        S = H{cUser.baseNode,iBand}(:,:,iUser) * V{iUser,1};

        for jUser = 1:SimParams.nUsers
            if ~isempty(V{jUser,1})
                if iUser ~= jUser
                    ifUser = SimStructs.userStruct{jUser,1};
                    I = H{ifUser.baseNode,iBand}(:,:,iUser) * V{jUser,1};
                    N = N + I * I';            
                end
            end
        end

        L = eye(size(N)) + (S * S') / N;
        SumCapacity(1,1) = SumCapacity(1,1) + log2(abs(det(L)));
        
    end
    
end

