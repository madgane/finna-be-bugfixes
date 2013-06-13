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
        
    case 'UMi'
        
        SimParams.sysConfig.ISD = 200;
        SimParams.sysConfig.carrierFreqGHz = 2.5;
        SimParams.sysConfig.BStransmitPwr_dBm = 41;

        SimParams.sysConfig.Kfactor.avg = 9;
        SimParams.sysConfig.Kfactor.std = 5;
        
        SimParams.sysConfig.shadowing.LOS = 3;
        SimParams.sysConfig.shadowing.NLOS = 4;
        SimParams.sysConfig.shadowing.OtoI = 7;
        
    case 'UMa'
        
        SimParams.sysConfig.ISD = 500;
        SimParams.sysConfig.carrierFreqGHz = 2.0;
        SimParams.sysConfig.BStransmitPwr_dBm = 46;
        
        SimParams.sysConfig.Kfactor.avg = 9;
        SimParams.sysConfig.Kfactor.std = 3.5;
        
        SimParams.sysConfig.shadowing.LOS = 4;
        SimParams.sysConfig.shadowing.NLOS = 6;
        SimParams.sysConfig.shadowing.OtoI = 0;
        
    case 'SMa'
        
        SimParams.sysConfig.ISD = 1299;
        SimParams.sysConfig.carrierFreqGHz = 2.0;
        SimParams.sysConfig.BStransmitPwr_dBm = 46;
        
        SimParams.sysConfig.Kfactor.avg = 9;
        SimParams.sysConfig.Kfactor.std = 7;

        SimParams.sysConfig.shadowing.LOS = 4;
        SimParams.sysConfig.shadowing.NLOS = 8;
        SimParams.sysConfig.shadowing.OtoI = 0;
                
    case 'RMa'
        
        SimParams.sysConfig.ISD = 1732;
        SimParams.sysConfig.carrierFreqGHz = 0.8;
        SimParams.sysConfig.BStransmitPwr_dBm = 46;
        
        SimParams.sysConfig.Kfactor.avg = 7;
        SimParams.sysConfig.Kfactor.std = 4;

        SimParams.sysConfig.shadowing.LOS = 4;
        SimParams.sysConfig.shadowing.NLOS = 8;
        SimParams.sysConfig.shadowing.OtoI = 0;
        
end

SimParams.sysConfig.systemBWHz = 10e6;SimParams.sysConfig.systemTones = 1024;
SimParams.sysConfig.NoisePwr_dBm = -174 + 10 * log10(SimParams.sysConfig.systemBWHz);

end

