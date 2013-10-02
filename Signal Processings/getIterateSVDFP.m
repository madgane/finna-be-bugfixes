
function [U1, D1, V1] = getIterateSVDFP(X,nIter,wordlength,fraclength)

X = sfi(X',wordlength,fraclength);

iIter = 1;
[nR,nC] = size(X.data);
Q = cell(nIter,1);R = cell(nIter,1);

while iIter <= nIter
    
    if iIter == 1
        [~,~,Q{iIter,1}, R{iIter,1}] = getGivensRotationFP(X,wordlength,fraclength);
    else
        [~,~,Q{iIter,1}, R{iIter,1}] = getGivensRotationFP(R{iIter - 1,1}',wordlength,fraclength);
    end
    iIter = iIter + 1;
end

U = sfi(eye(nC),wordlength,fraclength);
V = sfi(eye(nR),wordlength,fraclength);

for iIter = 1:nIter
    if mod(iIter - 1,2) == 0
        V = quantize(V * Q{iIter,1},1,wordlength,fraclength);
    else
        U = quantize(Q{iIter,1}' * U,1,wordlength,fraclength);
    end
end

if mod(nIter,2)
    D = R{nIter,1}';
else
    D = R{nIter,1};
end

V1 = V.data;U1 = U.data';D1 = D.data;
    
end
