function displayOutputs(SimParams, SimStructs)

switch SimParams.plotMode
    
    case 'SRA'
        
        SimResults.sumThrpt = sum(SimParams.Thrpt(:,:,end),2);
        SimResults.thrptFairness = sum(SimParams.fairness(:,:,end),2);
        SimParams.sumThrpt = SimResults.sumThrpt;
        
        figStruct.N = 1;figStruct.P = 'plot';
        figStruct.X = SimParams.snrIndex;figStruct.Y = SimParams.sumThrpt;
        
        plotFigure(figStruct);
        xlabel('SNR in dB');ylabel('sum rate in bits/sec/Hz');
        
        JainMean = mean(SimParams.Thrpt,2).^2;JainVar = var(SimParams.Thrpt,0,2);
        JainIndex_capacity = JainMean ./ (JainMean + JainVar);
        
        figStruct.N = 2;figStruct.P = 'plot';
        figStruct.X = SimParams.snrIndex;figStruct.Y = JainIndex_capacity;

%         plotFigure(figStruct);
%         xlabel('SNR in dB');ylabel('Rate Deviation across Users in bits/sec/Hz');
        
        JainMean = mean(SimParams.fairness,2).^2;JainVar = var(SimParams.fairness,0,2);
        JainIndex_utility = JainMean ./ (JainMean + JainVar);
        
        figStruct.N = 2;figStruct.P = 'plot';
        figStruct.X = SimParams.snrIndex;figStruct.Y = JainIndex_utility;

%         plotFigure(figStruct);
%         xlabel('SNR in dB');ylabel('Network Utility Deviation across Users');
        
        
    case 'QA'

        figStruct.N = 1;figStruct.P = 'plot';
        figStruct.X = 1:SimParams.nDrops;figStruct.Y = sum(squeeze(SimParams.QueueInfo.queueBackLogsOverTime(end,:,end,:)),1);
        
        plotFigure(figStruct);
        xlabel('Slot Index');ylabel('Queue Backlogs (pkts) over Time');grid on;
        
        %         plotFigure(1:SimParams.nDrops,std(squeeze(SimParams.QueueInfo.queueBackLogsOverTime(end,:,end,:)),1),5,'plot');
        %         xlabel('Slot Index');ylabel('{\sigma_Q} Queue Backlogs (pkts) over Time');grid on;
        
        %         plotFigure(SimParams.maxArrival,sum(squeeze(SimParams.QueueInfo.queueBackLogs(end,:,:)),1),6,'plot');
        %         xlabel('Average Arrival Rate');ylabel('Average Queue Size (pkts)');grid on;
        
    case 'STA'
        
        nT = 1e3;nPRB = 50;nREinPRB = 120;nTot = nT * nPRB * nREinPRB * 1e-6;

        figStruct.N = 1;figStruct.P = 'cdfplot';
        figStruct.Y = SimParams.Thrpt(1,:,1) * nTot;
        
        plotFigure(figStruct);
        xlabel('Throughput in Mbps');
        ylabel('CDF of Throughput in Mbps');
        
    case 'NRA'
        
        plotFigure(struct('Y',SimParams.sumRateInstant));
        
    case 'QInfo'
        
        clc;
        
        displaySystemDetails;
        displayChannel(SimParams,SimStructs);
        
        for iDrop = SimParams.nDrops:SimParams.nDrops
            displayQueues(SimParams,SimStructs,iDrop);
        end
        
    case 'QTimePlot'
        
        plotFigure(struct('Y',sum(squeeze(SimParams.QueueInfo.queueResiduesOverTime(end,:,end,:)))));
        
    case 'CPlot'
        
        displaySystemDetails;
        displayChannel(SimParams,SimStructs);
        displayQueues(SimParams,SimStructs);
        plotFigure(struct('Y',sum(cell2mat(SimParams.Debug.tempResource{3,1})),'N',1));
        xlabel('Iteration count');ylabel('Queue deviation in bits / channel use');
        
        plotFigure(struct('Y',sum(cell2mat(SimParams.Debug.tempResource{2,1})),'N',2));
        xlabel('Iteration count');ylabel('Sum rate in bits / channel use');        
        
    otherwise
        
        display('Simulation Completed without any display !');
        
end

end