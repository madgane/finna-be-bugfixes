
if SimParams.iSNR * SimParams.iDrop == 1    
    if strcmpi(SimParams.saveChannelInfo,'true')
        channelFile = sprintf('%s\\channelDumpFile.dat',SimParams.channelSaveFolder);
        fid = fopen(channelFile,'w');
        
        for sciUser = 1:SimParams.nUsers
            for sciBand = 1:SimParams.nBands
                for sciBase = 1:SimParams.nBases
                    for sciRow = 1:SimParams.nRxAntenna
                        for sciCol = 1:SimParams.nTxAntenna
                            fprintf(fid,'%f\t%f\t',real(SimStructs.linkChan{sciBase,sciBand}(sciRow,sciCol,sciUser)),imag(SimStructs.linkChan{sciBase,sciBand}(sciRow,sciCol,sciUser)));
                        end
                    end
                end
            end
        end
        
        fclose(fid);
    end    
end