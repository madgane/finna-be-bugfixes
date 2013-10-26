
function displayTimeBehaviour(matFileToLoad)

if ~nargin
    clc;clear all;
    globalCount = 1;
    load defaultTimeDomainFile.mat;
else
    load matFileToLoad;
end

legendString = cell(1,globalCount);

figLineWidth = {1};
figLineType = {'-','--','-.',':'};
figColor = {'b','g','r','m','c','k'};
figMarker = {'.','x','+','v','p','s','d','o'};

for iScheme = 1:globalCount
    
    SimParams = SimParamsCell{iScheme,1};
    SimStructs = SimStructsCell{iScheme,1};    
    displaySystemDetails;displayQueues(SimParams,SimStructs);
    
    fcIndex = mod(iScheme - 1,(length(figColor))) + 1;
    fmIndex = mod(iScheme - 1,(length(figMarker))) + 1;
    fltIndex = mod(iScheme - 1,(length(figLineType))) + 1;
    flwIndex = mod(iScheme - 1,(length(figLineWidth))) + 1;
    
    yValues = sum(squeeze(SimParams.QueueInfo.queueBacklogsOverTime(end,:,end,:)));
    plot(yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex});        
    
    hold all;
    legendString{1,iScheme} = SimParams.weightedSumRateMethod;

end

box on;
legend(legendString);

end
