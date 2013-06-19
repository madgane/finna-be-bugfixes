function [rssiMeasures, losMeasures] = getRSSIMeasures(SimParams,userPosition)

wrapModelLayout = 7;
nCites = getCellsOverLayout(SimParams.nTiers,1);
rssiMeasures = zeros(SimParams.nBases,wrapModelLayout);
losMeasures = cell(SimParams.nBases,wrapModelLayout);

for iWrapMode = 1:wrapModelLayout
    for iCite = 1:nCites
        for iSector = 1:SimParams.nSectors
            currentSite = (iCite - 1) * SimParams.nSectors + iSector;
            basePosition = SimParams.wrapCellLocArray(iCite,iWrapMode);
            
            separationM = abs(userPosition - basePosition);
                       
            [xRSSI, losMeasures{currentSite,iWrapMode}] = evaluateLTE_PL(SimParams,separationM,'false');
            antennaGain = getAntennaPatterGain(basePosition,userPosition,SimParams.sysConfig.layoutFeatures,iSector);
                       
            rssiMeasures(currentSite,iWrapMode) = xRSSI + antennaGain;
        end
    end
end

end