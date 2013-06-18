function [antPatternGain] = getAntennaPatterGain(varargin)

baseLoc = varargin{1,1};
userLoc = varargin{1,2};
layoutParams = varargin{1,3};

halfElAngle = 15;
halfAzAngle = 70;
minAntGain_dB = 20;

theta = angle(userLoc - baseLoc) * 180 / pi - layoutParams.layoutAngleFromEast;
phi = atan((layoutParams.hBS - layoutParams.hUT)/abs(userLoc - baseLoc)) * 180 / pi;

azAntGain = -min(12 * (theta / halfAzAngle)^2 , minAntGain_dB);
elAntGain = -min(12 * ((phi - layoutParams.antTilt) / halfElAngle)^2 , minAntGain_dB);

antPatternGain = -min(-(azAntGain + elAntGain),minAntGain_dB);
    
end



