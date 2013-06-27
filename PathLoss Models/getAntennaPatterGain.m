function [antPatternGain] = getAntennaPatterGain(varargin)

baseLoc = varargin{1,1};
userLoc = varargin{1,2};
layoutParams = varargin{1,3};
cSector = varargin{1,4};
nSectors = varargin{1,5};

limitAngle = 180;
halfElAngle = 15;
halfAzAngle = 70;
minAntGain_dB = 20;

cAngle = -90 + 120 * (cSector - 1);
theta = angle(userLoc - baseLoc) * 180 / pi - layoutParams.layoutAngleFromEast - cAngle;

if theta >= limitAngle
    theta = theta - limitAngle * 2;
end
if theta <= -limitAngle
    theta = theta + limitAngle * 2;
end

phi = atan((layoutParams.hBS - layoutParams.hBS)/abs(userLoc - baseLoc)) * 180 / pi;

switch nSectors
    case 1
        azAntGain = 0;
    case 3
        azAntGain = -min(12 * (theta / halfAzAngle)^2 , minAntGain_dB);
end

elAntGain = -min(12 * ((phi - layoutParams.antTilt) / halfElAngle)^2 , minAntGain_dB);

antPatternGain = -min(-(azAntGain + elAntGain),minAntGain_dB);

end



