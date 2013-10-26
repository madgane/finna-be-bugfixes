function displayOutputs(SimParams, SimStructs)

switch SimParams.plotMode
    
    case 'SRA'
        
        SimResults.sumThrpt = sum(SimParams.Thrpt(:,:,end),2);
        SimResults.thrptFairness = sum(SimParams.fairness(:,:,end),2);
        SimParams.sumThrpt = SimResults.sumThrpt;
        
        plotFigure(SimParams.snrIndex,SimParams.sumThrpt,1,'plot');
        xlabel('SNR in dB');ylabel('sum rate in bits/sec/Hz');
        
        JainMean = mean(SimParams.Thrpt,2).^2;JainVar = var(SimParams.Thrpt,0,2);
        JainIndex_capacity = JainMean ./ (JainMean + JainVar);
        
        %         plotFigure(SimParams.snrIndex,JainIndex_capacity,2,'plot');
        %         xlabel('SNR in dB');ylabel('Rate Deviation across Users in bits/sec/Hz');
        
        JainMean = mean(SimParams.fairness,2).^2;JainVar = var(SimParams.fairness,0,2);
        JainIndex_utility = JainMean ./ (JainMean + JainVar);
        
        %         plotFigure(SimParams.snrIndex,JainIndex_utility,3,'plot');
        %         xlabel('SNR in dB');ylabel('Network Utility Deviation across Users');
        
        
    case 'QA'
        
        plotFigure(1:SimParams.nDrops,sum(squeeze(SimParams.QueueInfo.queueBackLogsOverTime(end,:,end,:)),1),4,'plot');
        xlabel('Slot Index');ylabel('Queue Backlogs (pkts) over Time');grid on;
        
        %         plotFigure(1:SimParams.nDrops,std(squeeze(SimParams.QueueInfo.queueBackLogsOverTime(end,:,end,:)),1),5,'plot');
        %         xlabel('Slot Index');ylabel('{\sigma_Q} Queue Backlogs (pkts) over Time');grid on;
        
        %         plotFigure(SimParams.maxArrival,sum(squeeze(SimParams.QueueInfo.queueBackLogs(end,:,:)),1),6,'plot');
        %         xlabel('Average Arrival Rate');ylabel('Average Queue Size (pkts)');grid on;
        
    case 'STA'
        
        nT = 1e3;nPRB = 50;nREinPRB = 120;nTot = nT * nPRB * nREinPRB * 1e-6;
        
        hold all;
        plotFigure(SimParams.Thrpt(1,:,1) * nTot,1,1,'cdfplot');
        xlabel('Throughput in Mbps');
        ylabel('CDF of Throughput in Mbps');
        
    case 'NRA'
        
        plotFigure(1:SimParams.nDrops,SimParams.sumRateInstant,1,'plot');
        
    case 'QInfo'
        
        clc;
        
        displaySystemDetails;
        displayChannel(SimParams,SimStructs);
        
        for iDrop = SimParams.nDrops:SimParams.nDrops
            displayQueues(SimParams,SimStructs,iDrop);
        end
        
    otherwise
        
        display('Simulation Completed without any display !');
        
end

end