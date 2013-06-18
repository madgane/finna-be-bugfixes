function [SimParams, SimStructs] = configureLTEParams(SimParams,SimStructs)

plModel = char(SimParams.pathLossModel);
uscore_index = find(plModel == '_');
pathLossModel = plModel(uscore_index(1,1) + 1:end);

switch pathLossModel
    
    case 'InH'
        
        SimParams.sysConfig.ISD = 60;
        SimParams.sysConfig.carrierFreqGHz = 3.4;
        SimParams.sysConfig.BStransmitPwr_dBm = 21;
        
        SimParams.sysConfig.Kfactor.avg = 7;
        SimParams.sysConfig.Kfactor.std = 4;
                
        SimParams.sysConfig.shadowing.LOS = 3;
        SimParams.sysConfig.shadowing.NLOS = 4;
        SimParams.sysConfig.shadowing.OtoI = 0;
        
        SimParams.sysConfig.baseTerminalBG = 0;
        
        SimParams.sysConfig.layoutFeatures.minDistance = 3.0;
        SimParams.sysConfig.layoutFeatures.hUT = 2.5;
        SimParams.sysConfig.layoutFeatures.hBS = 6.0;
        
        
    case 'UMi'
        
        SimParams.sysConfig.ISD = 200;
        SimParams.sysConfig.carrierFreqGHz = 2.5;
        SimParams.sysConfig.BStransmitPwr_dBm = 41;

        SimParams.sysConfig.Kfactor.avg = 9;
        SimParams.sysConfig.Kfactor.std = 5;
        
        SimParams.sysConfig.shadowing.LOS = 3;
        SimParams.sysConfig.shadowing.NLOS = 4;
        SimParams.sysConfig.shadowing.OtoI = 7;
        
        SimParams.sysConfig.baseTerminalBG = 0;
        
        SimParams.sysConfig.layoutFeatures.minDistance = 10.0;
        SimParams.sysConfig.layoutFeatures.hUT = 2.5;
        SimParams.sysConfig.layoutFeatures.hBS = 10.0;

        
    case 'UMa'
        
        SimParams.sysConfig.ISD = 500;
        SimParams.sysConfig.carrierFreqGHz = 2.0;
        SimParams.sysConfig.BStransmitPwr_dBm = 46;
        
        SimParams.sysConfig.Kfactor.avg = 9;
        SimParams.sysConfig.Kfactor.std = 3.5;
        
        SimParams.sysConfig.shadowing.LOS = 4;
        SimParams.sysConfig.shadowing.NLOS = 6;
        SimParams.sysConfig.shadowing.OtoI = 17;
        
        SimParams.sysConfig.baseTerminalBG = 17;
        
        SimParams.sysConfig.layoutFeatures.minDistance = 25.0;
        SimParams.sysConfig.layoutFeatures.hUT = 1.5;
        SimParams.sysConfig.layoutFeatures.hBS = 25.0;

    case 'SMa'
        
        SimParams.sysConfig.ISD = 1299;
        SimParams.sysConfig.carrierFreqGHz = 2.0;
        SimParams.sysConfig.BStransmitPwr_dBm = 46;
        
        SimParams.sysConfig.Kfactor.avg = 9;
        SimParams.sysConfig.Kfactor.std = 7;

        SimParams.sysConfig.shadowing.LOS = 4;
        SimParams.sysConfig.shadowing.NLOS = 8;
        SimParams.sysConfig.shadowing.OtoI = 0;
        
        SimParams.sysConfig.baseTerminalBG = 17;

        SimParams.sysConfig.layoutFeatures.minDistance = 35.0;
        SimParams.sysConfig.layoutFeatures.hUT = 1.5;
        SimParams.sysConfig.layoutFeatures.hBS = 35.0;
        
    case 'RMa'
        
        SimParams.sysConfig.ISD = 1732;
        SimParams.sysConfig.carrierFreqGHz = 0.8;
        SimParams.sysConfig.BStransmitPwr_dBm = 46;
        
        SimParams.sysConfig.Kfactor.avg = 7;
        SimParams.sysConfig.Kfactor.std = 4;

        SimParams.sysConfig.shadowing.LOS = 4;
        SimParams.sysConfig.shadowing.NLOS = 8;
        SimParams.sysConfig.shadowing.OtoI = 0;
        
        SimParams.sysConfig.baseTerminalBG = 17;

        SimParams.sysConfig.layoutFeatures.minDistance = 35.0;
        SimParams.sysConfig.layoutFeatures.hUT = 1.5;
        SimParams.sysConfig.layoutFeatures.hBS = 35.0;

end

SimParams.sysConfig.cableLoss = 2;
SimParams.sysConfig.baseTerminalNF = 5;
SimParams.sysConfig.userTerminalNF = 7;
SimParams.sysConfig.userTerminalBG = 0;

SimParams.sysConfig.layoutFeatures.antTilt = 10;
SimParams.sysConfig.layoutFeatures.layoutAngleFromEast = 60.00;

SimParams.sysConfig.systemBWHz = 10e6;SimParams.sysConfig.systemTones = 1024;
SimParams.sysConfig.NoisePwr_dBm = -174 + 10 * log10(SimParams.sysConfig.systemBWHz);

end

