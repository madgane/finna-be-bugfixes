function [SimParams,SimStructs,varargout] = getPMatrix(SimParams,SimStructs,varargin)

if nargin == 2
    
    switch (SimParams.PrecodingMethod)
        
        case 'Best_ZF_Method'
            [SimParams,SimStructs] = getZFMatrix(SimParams,SimStructs);
            
        case 'Best_Leakage_Method'
            [SimParams,SimStructs] = getLeakageMatrix(SimParams,SimStructs);
            
        case 'Best_BF_Method'
            [SimParams,SimStructs] = getBFMatrix(SimParams,SimStructs);
            
        case 'Best_MZF_Method'
            [SimParams,SimStructs] = getMZFMatrix(SimParams,SimStructs);
            
        case 'Best_Network_Method'
            [SimParams,SimStructs] = getNetworkMatrix(SimParams,SimStructs);
            
        case 'Best_DP_Method'
            [SimParams,SimStructs] = getDPMatrix(SimParams,SimStructs);
            
        case 'Best_RobustMinSINR_Method'
            [SimParams,SimStructs] = getRobustBeamformers(SimParams,SimStructs);
            
        case 'Best_WSR_Method'
            [SimParams,SimStructs] = getWSRMatrix(SimParams,SimStructs);
            
        case 'Best_WMMSE_Method'
            [SimParams,SimStructs] = getWMMSEMatrix(SimParams,SimStructs);
            
        case 'Best_TDMZF_Method'
            [SimParams,SimStructs] = getTDMZFMatrix(SimParams,SimStructs);
            
        case 'Best_CZF_Method'
            [SimParams,SimStructs] = getCZFMatrix(SimParams,SimStructs);
            
        case 'Best_WSRM_SCP_Method'
            [SimParams,SimStructs] = getWSRM_SCP(SimParams,SimStructs);
            
        otherwise
            display('Unknown Precoder Design !');
            
    end
    
else
    
    switch (SimParams.PrecodingMethod)
        
        case 'Best_ZF_Method'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'ZF');
            
        case 'Best_MZF_Method'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'MZF');
            
        case 'Best_WMMSE_Method'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'WMMSE');
            
        case 'Best_TDMZF_Method'
                        
        case 'Best_CZF_Method'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'J-ZF');
            
        case 'Best_WSRM_SCP_Method'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'WSRM');
            
        otherwise
            display('Unknown Precoder Design !');
            
    end
    
end

