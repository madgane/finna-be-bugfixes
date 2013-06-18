
clc;
clear all;

nIters = 10000;
X = zeros(nIters,1);

for iIter = 1:nIters
    X(iIter,1) = getPointInRhombus(100,10);
end

% hold all;
% plot(X,'or');
% A = exp(sqrt(-1) * 2 * (pi / 3)) * X;
% plot(A,'ob');
% A = exp(sqrt(-1) * 4 * (pi / 3)) * X;
% plot(A,'om');
% 
layoutParams.hBS = 25;
layoutParams.hUT = 1.5;
layoutParams.antTilt = 15;
layoutParams.layoutAngleFromEast = 60;

Y = zeros(nIters,1);
for iIter = 1:nIters
    Y(iIter,1) = getAntennaPatterGain(0+sqrt(-1)*0,X(iIter,1),layoutParams);
end

hold all;
plot(angle(X) * 180 / pi,Y,'o')


