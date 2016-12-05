function b = isNoseBehindHead(el_X , el_Y , nose_x, nose_y,  cameraCenter)
    y1 = nose_y; x2 = cameraCenter(1);
    x1 = nose_x; y2 = cameraCenter(2);
    m = (y1 - y2)/(x1 - x2);
    k = (x1*y2 - x2*y1)/(x1 - x2);
    if nose_x < cameraCenter(1)
        xRange = nose_x:0.1:cameraCenter(1);
    else
        xRange = cameraCenter(1):0.1:nose_x;
    end
    yRange = xRange*m + k;
    L1 = [el_X; el_Y];
    L2 = [xRange;yRange];
    P = InterX(L1,L2);
    if ~isempty(P)
        P = find((abs(P(1,:)-nose_x))>2 | (abs(P(2,:)-nose_y))>2, 1);
    end
    b = ~isempty(P);
end