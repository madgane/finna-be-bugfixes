function [SimParams,SimStructs] = userLayoutGeneration(SimParams,SimStructs)

tierAngle = pi / 3;
hexSide = SimParams.sysConfig.ISD / (2 * sqrt(3));
nUsersPerCell = SimParams.nUsers / SimParams.nBases;

SimStructs.baseStruct{1,1}.loc = complex(0,0);
for iBase = 2:SimParams.nBases
    SimStructs.baseStruct{iBase,1}.loc = exp(sqrt(-1) * ((iBase - 1) * tierAngle)) * SimParams.sysConfig.ISD;
end

cUser = 0;
for iCell = 1:SimParams.nBases
    for iUser = 1:nUsersPerCell
        cUser = cUser + 1;
        currentLoc = complex(rand,rand) * hexSide;        
        SimStructs.userStruct{cUser,1}.loc = SimStructs.baseStruct{iCell,1}.loc + currentLoc;
    end
end
