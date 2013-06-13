
function [Q1, R1] = qr_householder(X)

[nRow,nCol] = size(X);

R = X;
Q = eye(nRow);

for iCol = 1:min(nCol,nRow)
    s = -sign(R(iCol,iCol));
    normX = norm(R(iCol:end,iCol));

    u1 = R(iCol,iCol) - s * normX;
    w = R(iCol:end,iCol)/(u1);
    w(1,1) = 1;
    tau = -s * u1 / normX;
    
    R(iCol:end,:) = R(iCol:end,:)-(tau*w)*(w'*R(iCol:end,:));
    Q(:,iCol:end) = Q(:,iCol:end)-(Q(:,iCol:end)*w)*(tau*w)';
end

Q1 = Q;R1 = R;

end


