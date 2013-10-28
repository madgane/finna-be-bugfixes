
function [varargout] = randomizeInitialSCApoint(varargin)

switch nargin
    case 2
        SimParams = varargin{1};
        SimStructs = varargin{2};
    case 3
        SimParams = varargin{1};
        SimStructs = varargin{2};
        currentBand = varargin{3};
    otherwise
        display('Unknown arguments !');
end

cH = SimStructs.linkChan;
nBases = SimParams.nBases;
nBands = SimParams.nBands;

usersPerCell = zeros(nBases,1);
cellUserIndices = cell(nBases,1);

for iBase = 1:nBases
    for iBand = 1:nBands
        cellUserIndices{iBase,1} = [cellUserIndices{iBase,1} ; SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}];
    end
    cellUserIndices{iBase,1} = unique(cellUserIndices{iBase,1});
    usersPerCell(iBase,1) = length(cellUserIndices{iBase,1});
end

nUsers = sum(usersPerCell);
maxRank = SimParams.maxRank;

W = cell(nUsers,nBands);
p = zeros(maxRank,nUsers,nBands);
q = zeros(maxRank,nUsers,nBands);
b = zeros(maxRank,nUsers,nBands);
M = ones(SimParams.nTxAntenna,maxRank,nUsers,nBands);

for iBand = 1:nBands
    for iBase = 1:nBases
        for iUser = 1:usersPerCell(iBase,1)
            cUser = cellUserIndices{iBase,1}(iUser,1);
            for iLayer = 1:maxRank
                R = SimParams.N * eye(SimParams.nRxAntenna);
                for jBase = 1:nBases
                    for jUser = 1:usersPerCell(jBase,1)
                        rUser = cellUserIndices{jBase,1}(jUser,1);
                        H = cH{jBase,iBand}(:,:,cUser);
                        R = R + H * M(:,:,rUser,iBand) * M(:,:,rUser,iBand)' * H';
                    end
                end
                H = cH{iBase,iBand}(:,:,cUser);
                W{cUser,iBand}(:,iLayer) = R \ (H * M(:,iLayer,cUser,iBand));
            end
        end
    end
end

for iBand = 1:nBands
    for iUser = 1:nUsers
        for iLayer = 1:maxRank
            
            cW = W{iUser,iBand}(:,iLayer);
            N = SimParams.N * norm(cW)^2;
            for jUser = 1:nUsers
                if jUser ~= iUser
                    for jLayer = 1:maxRank
                        N = N + norm(cW' * cH{SimStructs.userStruct{jUser,1}.baseNode,iBand}(:,:,iUser) * M(:,jLayer,jUser,iBand))^2;
                    end
                else
                    for jLayer = 1:maxRank
                        if iLayer ~= jLayer
                            N = N + norm(cW' * cH{SimStructs.userStruct{jUser,1}.baseNode,iBand}(:,:,iUser) * M(:,jLayer,jUser,iBand))^2;
                        end
                    end
                end
            end
            
            b(iLayer,iUser,iBand) = N;
            p(iLayer,iUser,iBand) = real(cW' * cH{SimStructs.userStruct{iUser,1}.baseNode,iBand}(:,:,iUser) * M(:,iLayer,iUser,iBand));
            q(iLayer,iUser,iBand) = imag(cW' * cH{SimStructs.userStruct{iUser,1}.baseNode,iBand}(:,:,iUser) * M(:,iLayer,iUser,iBand));
            
        end
        
    end
end

switch nargin
    case 2
        varargout{1} = p;varargout{2} = q;varargout{3} = b;
        if nargout == 4
            varargout{4} = W;
        end
    case 3
    varargout{1} = p(:,:,currentBand);varargout{2} = q(:,:,currentBand);varargout{3} = b(:,:,currentBand);
    if nargout == 4
        xW = cell(nUsers,1);
        for iUser = 1:nUsers
            xW{iUser,1} = W{iUser,currentBand};
        end
        varargout{4} = xW;
    end
end

end






