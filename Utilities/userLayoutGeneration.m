function [SimParams,SimStructs] = userLayoutGeneration(SimParams,SimStructs)

maxTiers = 3;
nSectors = 3;
tierAngle = pi / 3;
hexSide = SimParams.sysConfig.ISD / (2 * sqrt(3));
nUsersPerCell = SimParams.nUsers / SimParams.nBases;

totCells = 0;
for iTier = 1:maxTiers
    if iTier == 1
        nCells = 1;
    else
        nCells = iTier * (iTier + 1);
    end
    
    totCells = totCells + nCells;
    if SimParams.nBases <= totCells
        totTiers = iTier;
        break;
    end    
end

baseLocArray = zeros(totCells,1);

for iTier = 1:totTiers
    if iTier == 1
        nCells = 1;
    else
        nCells = iTier * (iTier + 1);
    end
    
    xCell = 0;
    for iCell = 1:nCells
        xCell = xCell + iCell;
        baseLocArray(iCell,1) = exp(sqrt(-1) * ((iCell - 1) * tierAngle)) * SimParams.sysConfig.ISD * (iTier - 1);
    end
end

for iBase = 1:SimParams.nBases
    SimStructs.baseStruct{iBase,1}.loc = baseLocArray(iBase,1);
    for iUser = 1:nUsersPerCell
        cUser = nUsersPerCell * (iBase - 1) + iUser;
        userLoc = getPointInRhombus(hexSide,SimParams.sysConfig.layoutFeatures.minDistance);
        userLoc = exp(sqrt(-1) * randi(nSectors,1,1) * 2 * pi / 3) * userLoc;
        SimStructs.userStruct{cUser,1}.loc = SimStructs.baseStruct{iBase,1}.loc + userLoc;
    end
end

SimParams.wrapExtCells = baseLocArray;
SimParams.userArrayLocs = zeros(SimParams.nUsers,1);
SimParams.baseArrayLocs = zeros(SimParams.nBases,1);

for iBase = 1:SimParams.nBases
    SimParams.baseArrayLocs(iBase,1) = SimStructs.baseStruct{iBase,1}.loc;
end

for iUser = 1:SimParams.nUsers
    SimParams.userArrayLocs(iUser,1) = SimStructs.userStruct{iUser,1}.loc;
end
