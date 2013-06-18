function [varargout] = evaluateLTE_PL(SimParams,separationM)

varargout = cell(1,nargout);
spLight = 3e8;isLOS = 'false';
plModel = char(SimParams.pathLossModel);
uscore_index = find(plModel == '_');
pathLossModel = plModel(uscore_index(1,1) + 1:end);

switch pathLossModel
    
    case 'InH'
        
        if separationM <= 18
            prob_LOS = 1;
        elseif separationM < 37
            prob_LOS = exp(-(separationM - 18)/27);
        else
            prob_LOS = 0.5;
        end
        
        if rand < prob_LOS
            isLOS = 'true';
            pathLoss_dB = 16.9 * log10(separationM) + 32.8 + 20 * log10(SimParams.sysConfig.carrierFreqGHz);
            shadowFading = randn * SimParams.sysConfig.shadowing.LOS;
        else
            isLOS = 'false';
            pathLoss_dB = 43.3 * log10(separationM) + 11.5 + 20 * log10(SimParams.sysConfig.carrierFreqGHz);
            shadowFading = randn * SimParams.sysConfig.shadowing.NLOS;
        end
        
    case 'UMi'
        
        
    case 'UMa'
        
        W = 20;h = 20;
        hUT = SimParams.sysConfig.layoutFeatures.hUT;
        hBS = SimParams.sysConfig.layoutFeatures.hBS;
        
        prob_LOS = min((18/separationM),1) * (1 - exp (-separationM / 63)) + exp (-separationM / 63);
        breakDistance = (4 * (hBS - 1) * (hUT - 1) * SimParams.sysConfig.carrierFreqGHz * 1e9) / spLight;
        
        if rand < prob_LOS
            isLOS = 'true';
            if separationM < breakDistance
                pathLoss_dB = 22.0 * log10(separationM) + 28.0 + 20 * log10(SimParams.sysConfig.carrierFreqGHz);
            else
                pathLoss_dB = 40 * log10(separationM) + 7.8 - 18.0 * log10(hBS - 1) - 18.0 * log10(hUT - 1) + 20 * log10(SimParams.sysConfig.carrierFreqGHz);
            end
            
            shadowFading = randn * SimParams.sysConfig.shadowing.LOS;
        else
            isLOS = 'false';
            pathLoss_dB = 161.04 - 7.1 * log10(W) + 7.5 * log10(h) - (24.37 - 3.7 * (h/hBS)^2) * log10(hBS) ...
                + (43.42 - 3.1 * log10(hBS)) * (log10(separationM) - 3) + 20 * log10(SimParams.sysConfig.carrierFreqGHz) - (3.2 * (log10(11.75 * hUT))^2 - 4.97);
            shadowFading = randn * SimParams.sysConfig.shadowing.NLOS;
        end
        
        vehiclePenetrationLoss = 9 + 5 * randn;
        shadowFading = shadowFading + vehiclePenetrationLoss;
        
    case 'SMa'
        
        
    case 'RMa'
        
        
end

shadowFading = shadowFading + randn * SimParams.sysConfig.shadowing.OtoI;
powerCompensation = SimParams.sysConfig.BStransmitPwr_dBm + SimParams.sysConfig.userTerminalBG + SimParams.sysConfig.baseTerminalBG;
otherNoise = SimParams.sysConfig.baseTerminalNF + SimParams.sysConfig.userTerminalNF + SimParams.sysConfig.NoisePwr_dBm + SimParams.sysConfig.cableLoss;
varargout{1,1} = powerCompensation - pathLoss_dB - shadowFading - 10 * log10(SimParams.sysConfig.systemTones) - otherNoise;

if nargout > 1
    varargout{1,2} = isLOS;
end
