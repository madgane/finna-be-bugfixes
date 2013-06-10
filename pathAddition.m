
if isunix
    cd '~';cd 'furry-octo-hipster';
    addpath('..//simResults//');
else
    addpath('..\\simResults\\');
end

processID = feature('GetPid');
display(processID);

if isempty(strfind(path,'cvx'))
    
    currentDIR = dir('.');
    currentString = sprintf('%s',pwd);
    addpath(currentString);

    for iDirectory = 1:length(currentDIR)
       if currentDIR(iDirectory,1).isdir

           contWithDirectory = 1;
           if strcmp(currentDIR(iDirectory,1).name(1,1),'.')
               contWithDirectory = 0;
           end

           if contWithDirectory
               if isunix
                   currentString = sprintf('%s//%s',pwd,currentDIR(iDirectory,1).name);
               else
                   currentString = sprintf('%s\\%s',pwd,currentDIR(iDirectory,1).name);
               end
               addpath(currentString);
           end

       end
    end

    cvx_setup;
end

cvx_quiet('false');
cvx_solver('sedumi');

if isunix
    SimParams.outFile = sprintf('..//simResults//%s',SimParams.outFile);
else
    SimParams.outFile = sprintf('..\\simResults\\%s',SimParams.outFile);
end
