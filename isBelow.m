function bool_res = isBelow(mCoeff, kCoeff, x, y)
    % y = (x*y1 - x*y2 + x1*y2 - x2*y1)/(x1 - x2);
    % y = x*(y1 - y2)/(x1 - x2) + (x1*y2 - x2*y1)/(x1 - x2);
    yy = x*mCoeff+kCoeff;
    if(yy>y)
       bool_res = true;
    else
       bool_res = false;
    end
end

