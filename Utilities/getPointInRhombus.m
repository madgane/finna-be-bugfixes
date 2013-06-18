
function userLoc = getPointInRhombus(hexSide,minDistance)

while 1
    userLoc = complex(rand,rand) * hexSide;
    if abs(userLoc) > minDistance
        break;
    end        
end

K = tan(pi/6);K1 = sec(pi/6);
skewMat = [K1 -K ; 0 1];
userLoc = skewMat * [real(userLoc) ; imag(userLoc)];
userLoc = userLoc(1,1) + sqrt(-1) * userLoc(2,1);

end