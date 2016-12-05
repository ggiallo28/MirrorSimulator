function [ P1, P2 ] = getTangentLine(el_X,el_Y, cameraCenter)
% y = (x*y1 - x*y2 + x1*y2 - x2*y1)/(x1 - x2);
% y = x*(y1 - y2)/(x1 - x2) + (x1*y2 - x2*y1)/(x1 - x2);
% alpha = (y1 - y2)/(x1 - x2)
% c = (x1*y2 - x2*y1)/(x1 - x2)
alpha = zeros(1,size(el_X,2));
for i=1:size(el_X,2)
    y1 = cameraCenter(2);
    y2 = el_Y(i);
    x1 = cameraCenter(1);
    x2 = el_X(i);
    alpha(i) = (y1 - y2)/(x1 - x2);
end
P1(1) = el_X(alpha == min(alpha));
P1(2) = el_Y(alpha == min(alpha));
P2(1) = el_X(alpha == max(alpha));
P2(2) = el_Y(alpha == max(alpha));
end

