
Queues = zeros(SimParams.nUsers,1);
txPkts = zeros(SimParams.nUsers,SimParams.nBands);

for iUser = 1:SimParams.nUsers
    Queues(iUser,1) = SimStructs.userStruct{iUser,1}.trafficHistory.pktArrival(1,1);
    txPkts(iUser,:) = squeeze(SimParams.Debug.resAllocation(end,:,iUser,end));
end

servedPkts = sum(txPkts,2);
Qdeviation = sum(max((Queues - servedPkts),0));

display(Qdeviation);
[txPkts Queues servedPkts]

