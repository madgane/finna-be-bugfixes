
fprintf('\n');
display('Displaying User Queue status');
display('----------------------------');

printLatexScript = 'false';
Queues = zeros(SimParams.nUsers,1);
txPkts = zeros(SimParams.nUsers,SimParams.nBands);

for iUser = 1:SimParams.nUsers
    Queues(iUser,1) = SimStructs.userStruct{iUser,1}.trafficHistory.pktArrival(1,1);
    txPkts(iUser,:) = squeeze(SimParams.Debug.resAllocation(end,:,iUser,end));
end

servedPkts = sum(txPkts,2);
Qdeviation = sum(max((Queues - servedPkts),0));

QueueMatrix = [txPkts Queues servedPkts];

display(QueueMatrix);
fprintf('Queue Deviation - %f \n',Qdeviation);

% Displaying in Latex import format

if strcmp(printLatexScript,'true')
    
    [nRows,nCols] = size(QueueMatrix);

    for iRow = 1:nRows
        for iCol = 1:nCols
            fprintf('& %3.2f \t',QueueMatrix(iRow,iCol));
        end
        fprintf('\n');
    end
    
end

