function [nose_x, nose_y, rot] = rotateNose(nose, angleLeftMirror, angleRightMirror, mirrorsCenter, N, index )
    anglePoint = rad2deg(atan2((nose(2)-mirrorsCenter(2)),(nose(1)-mirrorsCenter(1))));
    if anglePoint<0
        anglePoint = -anglePoint;
    end
    v = genMirroring(angleLeftMirror, angleRightMirror, anglePoint, N);
    [nose_x, nose_y] = rotate(nose(1), nose(2) ,mirrorsCenter(1), mirrorsCenter(2), v(1,index), 'z'); 
    rot = v(1,index);
end

