function [totalRateAttained,QDeviation,rateVector] = performMockReception(SimParams,SimStructs,V,U,W,iBand)

% for iUser = 1:SimParams.nUsers
%     
%     if ~isempty(V{iUser,1})
%         
%         cUser = SimStructs.userStruct{iUser,1};    
%         S = U{iUser,1}' * H{cUser.baseNode,iBand}(:,:,iUser) * V{iUser,1};
%         N = SimParams.N * eye(SimParams.nRxAntenna) * (U{iUser,1}' * U{iUser,1});
% 
%         for jUser = 1:SimParams.nUsers
%             if ~isempty(V{jUser,1})
%                 if iUser ~= jUser
%                     ifUser = SimStructs.userStruct{jUser,1};
%                     I = U{iUser,1}' * H{ifUser.baseNode,iBand}(:,:,iUser) * V{jUser,1};
%                     N = N + I * I';            
%                 end
%             end
%         end
% 
%         L = eye(SimParams.nRxAntenna) + N \ (S * S');
%         SumCapacity = SumCapacity + real(log2(det(L)));
%         if nargout > 1 
%             rateVector(iUser,1) = real(log2(det(L)));
%             QDeviation(iUser,1) = max((SimStructs.userStruct{iUser,1}.weighingFactor - rateVector(iUser,1)),0);            
%         end        
% 
%     end
%     
% end

xParams = SimParams;
xStructs = SimStructs;
QueuedPkts = zeros(xParams.nUsers,1);
nStreams = min(SimParams.maxRank,SimParams.nRxAntenna);

for iBase = 1:xParams.nBases
    cBase = xStructs.baseStruct{iBase,1};
    
    cBase.P{iBand,1} = [];
    cBase.assignedUsers{iBand,1} = [];
    cBase.assignedStreams{iBand,1} = [];
    
    for iUser = 1:length(cBase.linkedUsers)
        cUserIndex = cBase.linkedUsers(iUser,1);
        cUser = xStructs.userStruct{cUserIndex,1};
        cBase.P{iBand,1} = [cBase.P{iBand,1} , V{cUserIndex,1}];
        cUser.W{iBand,1} = U{cUserIndex,1};
        
        xStreams = (1:nStreams)';
        cBase.assignedUsers{iBand,1} = [cBase.assignedUsers{iBand,1} ; repmat(cUserIndex,length(xStreams),1)];
        cBase.assignedStreams{iBand,1} = [cBase.assignedStreams{iBand,1} ; xStreams];
        xStructs.userStruct{cUserIndex,1} = cUser;
    end
    
    xStructs.baseStruct{iBase,1} = cBase;
end

for iUser = 1:xParams.nUsers
    QueuedPkts(iUser,1) = xStructs.userStruct{iUser,1}.weighingFactor;
end

[xParams,~] = performDummyReception(xParams,xStructs,iBand);

rateVector = xParams.Debug.privateExchanges.resAllocation(iBand,:)';
QDeviation = max(QueuedPkts - xParams.Debug.privateExchanges.resAllocation(iBand,:)',0);
totalRateAttained = sum(rateVector);


