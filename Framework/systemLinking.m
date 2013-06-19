function [SimParams,SimStructs] = systemLinking(SimParams,SimStructs)
    
for iBase = 1:SimParams.nBases
    SimStructs.baseStruct{iBase,1}.linkedUsers = [];
end

xBases = 1:SimParams.nBases;
[~,maxI] = max(SimParams.PL_Profile,[],1);

for iUser = 1:SimParams.nUsers
    cNode = maxI(1,iUser);
    SimStructs.userStruct{iUser,1}.baseNode = cNode;
    SimStructs.userStruct{iUser,1}.neighNode = find(xBases ~= cNode);
    SimStructs.baseStruct{cNode,1}.linkedUsers = [SimStructs.baseStruct{cNode,1}.linkedUsers ; iUser];
end

end