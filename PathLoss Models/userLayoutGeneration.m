function [SimParams,SimStructs] = userLayoutGeneration(SimParams,SimStructs)

debugRSSI = zeros(SimParams.nUsers,1);
userLocIndices = cell(SimParams.nBases,1);
nCites = getCellsOverLayout(SimParams.nTiers,1);
nUsersOverCell = floor(SimParams.nUsers / SimParams.nBases);

xUser = 0;
if SimParams.nSectors == 3
    hexSide = SimParams.sysConfig.ISD / 3;
else
    hexSide = SimParams.sysConfig.ISD / 2;
end
eastRotRad = SimParams.sysConfig.layoutFeatures.layoutAngleFromEast * pi / 180;

for iCite = 1:nCites
    citeLocation = SimParams.wrapCellLocArray(iCite,1);
    for iSector = 1:SimParams.nSectors
        
        cCite = (iCite - 1) * SimParams.nSectors + iSector;
        while (length(userLocIndices{cCite,1}) < nUsersOverCell)
            
            reIterate = 1;
            userPosition = 0;
            while reIterate
                switch SimParams.nSectors
                    case 1
                        userPosition = getPointInRhombus(hexSide,iSector,eastRotRad,0);
                    case 3
                        userPosition = getPointInRhombus(hexSide,iSector,eastRotRad,1);
                end
                
                if abs(userPosition) > SimParams.sysConfig.layoutFeatures.minDistance
                    if abs(userPosition) > SimParams.sysConfig.layoutFeatures.maxDistance
                        reIterate = 1;
                    else
                        reIterate = 0;
                    end
                else
                    reIterate = 1;
                end
                
            end
            
            userPosition = citeLocation + userPosition;
            [rssiMeasures, losMeasures] = getRSSIMeasures(SimParams,userPosition);
            
            [maxRSSI, maxRSSIindex] = max(rssiMeasures,[],2);
            [sortV,sortI] = sort(maxRSSI,'descend');
            
            if length(userLocIndices{sortI(1,1),1}) < nUsersOverCell
                xUser = xUser + 1;
                userLocIndices{sortI(1,1),1} = [userLocIndices{sortI(1,1),1} userPosition];
                SimStructs.userStruct{xUser,1}.phyParams.location = userPosition;
                
                SimStructs.userStruct{xUser,1}.phyParams.listedCites = sortI(1:(SimParams.nNeighbors + 1),1);
                SimStructs.userStruct{xUser,1}.phyParams.listedRSSI = sortV(1:(SimParams.nNeighbors + 1),1);
                
                for iNeighbor = 1:(SimParams.nNeighbors + 1)
                    SimStructs.userStruct{xUser,1}.losFading{sortI(iNeighbor,1),1} = losMeasures{sortI(iNeighbor,1),maxRSSIindex(sortI(iNeighbor,1))};
                end
                
                SimStructs.userStruct{xUser,1}.phyParams.restOfIF = 10 * log10(sum(10.^(sortV((SimParams.nNeighbors + 2):end,1)./10)));
                
                linRSSI = 10.^(sortV./10);
                debugRSSI(xUser,1) = 10 * log10(linRSSI(1,1) / (SimParams.N + sum(linRSSI(2:end,1))));
                
            end
        end
    end
end

end

