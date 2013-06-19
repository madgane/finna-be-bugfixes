function [SimParams,SimStructs] = userLayoutGeneration(SimParams,SimStructs)

debugRSSI = zeros(SimParams.nUsers,1);
userLocIndices = cell(SimParams.nBases,1);
nCites = getCellsOverLayout(SimParams.nTiers,1);
nUsersOverCell = floor(SimParams.nUsers / SimParams.nBases);

xUser = 0;
hexSide = SimParams.sysConfig.ISD / 3;
eastRotRad = SimParams.sysConfig.layoutFeatures.layoutAngleFromEast * pi / 180;

for iCite = 1:nCites
    citeLocation = SimParams.wrapCellLocArray(iCite,1);
    for iSector = 1:SimParams.nSectors
        
        cCite = (iCite - 1) * SimParams.nSectors + iSector;
        while (length(userLocIndices{cCite,1}) < nUsersOverCell)
            
            userPosition = 0;
            while abs(userPosition) < SimParams.sysConfig.layoutFeatures.minDistance
                userPosition = getPointInRhombus(hexSide,iSector,eastRotRad,1);
            end
            
            userPosition = citeLocation + userPosition;
            [rssiMeasures, losMeasures] = getRSSIMeasures(SimParams,userPosition);
            
            [maxRSSI, maxRSSIindex] = max(rssiMeasures,[],2);
            [sortV,sortI] = sort(maxRSSI,'descend');
            
            if length(userLocIndices{sortI(1,1),1}) < nUsersOverCell
                xUser = xUser + 1;
                userLocIndices{sortI(1,1),1} = [userLocIndices{sortI(1,1),1} userPosition];
                SimStructs.userStruct{xUser,1}.phyParams.location = userPosition;
                SimStructs.userStruct{xUser,1}.phyParams.servingBS = sortI(1,1);
                SimStructs.userStruct{xUser,1}.phyParams.servingRSSI = sortV(1,1);
                
                SimStructs.userStruct{xUser,1}.phyParams.neighborCites = sortI(2:(SimParams.nNeighbors + 1),1);
                SimStructs.userStruct{xUser,1}.phyParams.neighborRSSI = sortV(2:(SimParams.nNeighbors + 1),1);
                SimStructs.userStruct{xUser,1}.losFading{sortI(1,1),1} = losMeasures{sortI(1,1),maxRSSIindex(sortI(1,1))};
                
                for iNeighbor = 1:SimParams.nNeighbors
                    SimStructs.userStruct{xUser,1}.losFading{sortI((1 + iNeighbor),1),1} = losMeasures{sortI((1 + iNeighbor),1),maxRSSIindex(sortI((1 + iNeighbor),1))};
                end
                
                SimStructs.userStruct{xUser,1}.phyParams.restOfIF = 10 * log10(sum(10.^(sortV((SimParams.nNeighbors + 2):end,1)./10)));
                
                linRSSI = 10.^(sortV./10);
                debugRSSI(xUser,1) = 10 * log10(linRSSI(1,1) / (SimParams.N + sum(linRSSI(2:end,1))));
                
            end
        end
    end
end

end

% hexSide = SimParams.sysConfig.ISD * SimParams.nTiers;
% 
% for iUser = 1:SimParams.nUsers
%     
%     userPosition = getPointInRhombus(hexSide,0,eastRotRad,0);
%     [rssiMeasures, losMeasures] = getRSSIMeasures(SimParams,userPosition);
%     
%     [maxRSSI, maxRSSIindex] = max(rssiMeasures,[],2);
%     [sortV,sortI] = sort(maxRSSI,'descend');
%     
%     if length(userLocIndices{sortI(1,1),1}) < nUsersOverCell
%         xUser = xUser + 1;
%         userLocIndices{sortI(1,1),1} = [userLocIndices{sortI(1,1),1} userPosition];
%         SimStructs.userStruct{xUser,1}.phyParams.location = userPosition;
%         SimStructs.userStruct{xUser,1}.phyParams.servingBS = sortI(1,1);
%         SimStructs.userStruct{xUser,1}.phyParams.servingRSSI = sortV(1,1);
%         
%         SimStructs.userStruct{xUser,1}.phyParams.neighborCites = sortI(2:(SimParams.nNeighbors + 1),1);
%         SimStructs.userStruct{xUser,1}.phyParams.neighborRSSI = sortV(2:(SimParams.nNeighbors + 1),1);
%         SimStructs.userStruct{xUser,1}.losFading{sortI(1,1),1} = losMeasures{sortI(1,1),maxRSSIindex(sortI(1,1))};
%         
%         for iNeighbor = 1:SimParams.nNeighbors
%             SimStructs.userStruct{xUser,1}.losFading{sortI((1 + iNeighbor),1),1} = losMeasures{sortI((1 + iNeighbor),1),maxRSSIindex(sortI((1 + iNeighbor),1))};
%         end
%         
%         SimStructs.userStruct{xUser,1}.phyParams.restOfIF = 10 * log10(sum(10.^(sortV((SimParams.nNeighbors + 2):end,1)./10)));
%         
%         linRSSI = 10.^(sortV./10);
%         debugRSSI(xUser,1) = 10 * log10(linRSSI(1,1) / (SimParams.N + sum(linRSSI(2:end,1))));
%         
%     end
%     
%     
% end
% 
% end

% hold all;
% plot(SimParams.wrapCellLocArray(:,1),'o');
% for iCell = 1:length(SimParams.wrapCellLocArray(:,1))
%     plot(userLocIndices{iCell,1},'.');
% end
% 
% cdfplot(debugRSSI);

