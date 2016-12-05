clc; clear all; close all;
%% This script can be used to know how many faces could be visibile by changing the parameters of the system
%% Costants: These parameters are used to describe the parts of the system, are fixed and they refers to the dimensions of the mirror system in TU Wien.
lengthMirrors = 65;         %cm
mirror2Pivot = 2.5;         %cm Ci sta una struttura ad angolo retto che collega lo specchio all'ase di rotazione, abbiamo lungo l'asse dello specchio uno sfasamento di 2.5 cm
offset2Mirror = 3;          %cm ed uno sfasamento di 3 cm in avanti.
pivot2pivot = 28.5;         %cm Distanza tra i due perni
% Those calculation are useful because of the modeling of the system. The problem is that there is an offset between the pivot and the starting point of the reflective surface.
% Here we take in account of this offset. The distance between the pivot and the starting of reflective surface create an right-angled triangle. So first we calculate the Hypotenuse 
% of this triangle, and than the inner angles. 
lengthPivotMirrors = lengthMirrors+mirror2Pivot;
lengthHypotenuse = sqrt(offset2Mirror^2 + mirror2Pivot^2);
angle2Pivot = rad2deg(asin(offset2Mirror/lengthHypotenuse));
%% Parameters: You can change those parameters in order to check a different configuration.
angleRightMirror = 45;                       %deg
angleLeftMirror = 180-angleRightMirror;      %deg
distanceCamera = 170;                        %cm
fovHCamera = 50;                             %deg
isIdealCamera = true;                        %the value fovHCamera is used only if the value isIdealCamera is true, otherwise we consider a camera with variable field of view
%% HEAD MODEL
lengthHead = 15;            %cm % vertical 2*radius
widthHead = 10;             %cm % horizontal 2*radius
headPosx = 0.5*pivot2pivot; %cm
headPosy = 14;              %cm
%% Logic
angleBetweenMirror = abs(angleRightMirror-angleLeftMirror);  % Angle betwwen the two mirrors
N = 360/angleBetweenMirror; % FORMULA 1: Number of reflections produced, those visible will be a subset. To use as verification.
% Position of the inner extremes of the two mirrors with respect to the axis center placed on the left pivot. 
inPointLeft(1) =  lengthHypotenuse*cos(deg2rad(angleLeftMirror-angle2Pivot));                   %X 
inPointLeft(2) =  lengthHypotenuse*sin(deg2rad(angleLeftMirror-angle2Pivot));                   %Y
inPointRight(1) = pivot2pivot + lengthHypotenuse*cos(deg2rad(angleRightMirror+angle2Pivot));    %X
inPointRight(2) = lengthHypotenuse*sin(deg2rad(angleRightMirror+angle2Pivot));                  %Y
distanceBetweenMirror(1) = pdist([inPointLeft;inPointRight],'euclidean'); %shortest
% Position of the externals extremes
exPointLeft(1) = inPointLeft(1) + lengthMirrors*cos(deg2rad(angleLeftMirror));                  %X
exPointLeft(2) = inPointLeft(2) + lengthMirrors*sin(deg2rad(angleLeftMirror));                  %Y
exPointRight(1) = inPointRight(1) + lengthMirrors*cos(deg2rad(angleRightMirror));               %X
exPointRight(2) = inPointRight(2) + lengthMirrors*sin(deg2rad(angleRightMirror));               %Y
distanceBetweenMirror(2) = pdist([exPointLeft;exPointRight],'euclidean');
% Slope of the two lines that go through the two mirrors, we need this two lines in order to obtain the position of the rotational center.
mLeft = ((-inPointLeft(2)) - (-exPointLeft(2)))/(inPointLeft(1) - exPointLeft(1));
kLeft = (inPointLeft(1)*(-exPointLeft(2)) - exPointLeft(1)*(-inPointLeft(2)))/(inPointLeft(1) - exPointLeft(1));
mRight = ((-inPointRight(2)) - (-exPointRight(2)))/(inPointRight(1) - exPointRight(1));
kRight = (inPointRight(1)*(-exPointRight(2)) - exPointRight(1)*(-inPointRight(2)))/(inPointRight(1) - exPointRight(1));
[mirrorsCenter(1),mirrorsCenter(2)] = inters2rette(mLeft,kLeft,mRight,kRight);
%% CAMERA POSITION
% If we have and ideal camera me consider the distanceCamera as parameter, otherwise we adjust the position of the camera according to the current
% angle betwwen the mirrors
cameraCenter(1) = 0.5*pivot2pivot; 
if isIdealCamera
    cameraCenter(2) = -distanceCamera; 
else
    cameraCenter(2) = -(exPointRight(2) +0.5*distanceBetweenMirror(2)*sqrt(1/cos(deg2rad(180-90-0.5*fovHCamera))^2-1));
end
% Uncomment to align the head and the camera with the mirror rotation center
% headPosx = mirrorsCenter(1);
% cameraCenter(1) = mirrorsCenter(1);

% Angle of the line through the center of the mirrors and the center of the head, is used to figure out how to rotate
angleHead = rad2deg(atan2((-headPosy-mirrorsCenter(2)),(headPosx-mirrorsCenter(1))));
if angleHead<0
    angleHead = -angleHead;
end
vH = genMirroring(angleLeftMirror, angleRightMirror, angleHead, N);

[el_x, el_y, nose] = genHead(headPosx, -headPosy, widthHead, lengthHead);
el_X = zeros(floor(N)-1,size(el_x,2));
el_Y = zeros(floor(N)-1,size(el_y,2));
distOR = sqrt((mirrorsCenter(1)-mean(el_x))^2 + (mirrorsCenter(2)-mean(el_y))^2);
% Compute Rotations
for i=1:size(el_x,2)
    anglePoint = rad2deg(atan2((el_y(i)-mirrorsCenter(2)),(el_x(i)-mirrorsCenter(1))));
    if anglePoint<0
        anglePoint = -anglePoint;
    end
    v = genMirroring(angleLeftMirror, angleRightMirror, anglePoint, N);
    for j=1:floor(N)-1
       [el_X(j,i),el_Y(j,i)] = rotate(el_x(i),el_y(i),mirrorsCenter(1),mirrorsCenter(2),v(1,j),'z'); 
    end
end
% Check Rotations
for j=1:floor(N)-1
    new_center_x = mean(el_X(j,:));
    new_center_y = mean(el_Y(j,:));
    distROT = sqrt((mirrorsCenter(1)-mean(el_X(j,:)))^2 + (mirrorsCenter(2)-mean(el_Y(j,:)))^2);
    assert(abs(distOR-distROT)<exp(-10));
end
%% Plot we have to use -Y. plot([x1 x2], [y1 y2])
figure,hold on,axis([-distanceCamera distanceCamera -distanceCamera distanceCamera])
plot([0 pivot2pivot], [0 0],'b');
plot([0 inPointLeft(1)], [0 -inPointLeft(2)],'b');
plot([pivot2pivot inPointRight(1)], [0 -inPointRight(2)],'b');
plot([inPointRight(1) exPointRight(1)], [-inPointRight(2) -exPointRight(2)],'b');
plot([inPointLeft(1) inPointRight(1)], [-inPointLeft(2) -inPointRight(2)],'b');
plot([inPointLeft(1) exPointLeft(1)], [-inPointLeft(2) -exPointLeft(2)],'b');

scatter(mirrorsCenter(1),mirrorsCenter(2));
scatter(cameraCenter(1),cameraCenter(2));
plot(el_x,el_y,'b')
plot([min(el_x) cameraCenter(1)], [median(el_y) cameraCenter(2)],'b');
plot([max(el_x) cameraCenter(1)], [median(el_y) cameraCenter(2)],'b');

for i=1:size(v,2);
    plot(el_X(i,:),el_Y(i,:));
    [p1, p2] = getTangentLine(el_X(i,:),el_Y(i,:), cameraCenter);
    plot([p1(1) cameraCenter(1)], [p1(2) cameraCenter(2)]);
    plot([p2(1) cameraCenter(1)], [p2(2) cameraCenter(2)]);
end
%% Drawing the scene
% y = (x*y1 - x*y2 + x1*y2 - x2*y1)/(x1 - x2);
% y = x*(y1 - y2)/(x1 - x2) + (x1*y2 - x2*y1)/(x1 - x2);
y1 = -inPointLeft(2); y2 = -inPointRight(2);
x1 = inPointLeft(1); x2 = inPointRight(1);
mCoeff = (y1 - y2)/(x1 - x2);
kCoeff = (x1*y2 - x2*y1)/(x1 - x2);
xRange = x1:0.1:x2;
yRange = xRange*mCoeff + kCoeff;

XRange = zeros(floor(N)-1,size(xRange,2));
YRange = zeros(floor(N)-1,size(yRange,2));
distOR = sqrt((mirrorsCenter(1)-mean(xRange))^2 + (mirrorsCenter(2)-mean(yRange))^2);
for i=1:size(xRange,2)
    anglePoint = rad2deg(atan2((yRange(i)-mirrorsCenter(2)),(xRange(i)-mirrorsCenter(1))));
    if anglePoint<0
        anglePoint = -anglePoint;
    end
    v = genMirroring(angleLeftMirror, angleRightMirror, anglePoint, N);
    for j=1:floor(N)-1
       [XRange(j,i),YRange(j,i)] = rotate(xRange(i),yRange(i),mirrorsCenter(1),mirrorsCenter(2),v(1,j),'z');
    end 
end

for j=1:floor(N)-1
    plot(XRange(j,:),YRange(j,:))
    new_center_x = mean(XRange(j,:));
    new_center_y = mean(YRange(j,:));
    distROT = sqrt((mirrorsCenter(1)-mean(XRange(j,:)))^2 + (mirrorsCenter(2)-mean(YRange(j,:)))^2);
    assert(abs(distOR-distROT)<exp(-10));
end
hold off;
%% Plot full setup
fig = figure; hold on,axis([-distanceCamera distanceCamera -distanceCamera distanceCamera])
plot([0 pivot2pivot], [0 0],'b');
plot([0 inPointLeft(1)], [0 -inPointLeft(2)],'b');
plot([pivot2pivot inPointRight(1)], [0 -inPointRight(2)],'b');
plot([inPointRight(1) exPointRight(1)], [-inPointRight(2) -exPointRight(2)],'b');
plot([inPointLeft(1) inPointRight(1)], [-inPointLeft(2) -inPointRight(2)],'b');
plot([inPointLeft(1) exPointLeft(1)], [-inPointLeft(2) -exPointLeft(2)],'b');

scatter(mirrorsCenter(1),mirrorsCenter(2));
scatter(cameraCenter(1),cameraCenter(2));
plot(el_x,el_y,'b')
plot([min(el_x) cameraCenter(1)], [median(el_y) cameraCenter(2)],'b');
plot([max(el_x) cameraCenter(1)], [median(el_y) cameraCenter(2)],'b');

% Line through the center of the camera and the outer end of the left mirror
y1 = -exPointLeft(2); x2 = cameraCenter(1);
x1 = exPointLeft(1); y2 = cameraCenter(2);
mCoeff2Left = (y1 - y2)/(x1 - x2);
kCoeff2Left = (x1*y2 - x2*y1)/(x1 - x2);
xRange2Left = -200:0.1:200;
yRange2Left = xRange2Left*mCoeff2Left + kCoeff2Left;
if ~isIdealCamera
    rad_fovHCamera = deg2rad(fovHCamera);
    mCoeff2Left_tmp = tan(pi/2-0.5*rad_fovHCamera)*sign(mCoeff2Left);
    if abs(mCoeff2Left_tmp) > abs(mCoeff2Left)
        mCoeff2Left = mCoeff2Left_tmp;
        yRange2Left = mCoeff2Left*(xRange2Left - cameraCenter(1)) + cameraCenter(2);  
        kCoeff2Left = (xRange2Left(1)*yRange2Left(2) - xRange2Left(2)*yRange2Left(1))/(xRange2Left(1) - xRange2Left(2));
    end
end
plot(xRange2Left,yRange2Left,'k');

% Line through the center of the camera and the outer end of the right mirror
y1 = -exPointRight(2); x2 = cameraCenter(1);
x1 = exPointRight(1); y2 = cameraCenter(2);
mCoeff2Right = (y1 - y2)/(x1 - x2);
kCoeff2Right = (x1*y2 - x2*y1)/(x1 - x2);
xRange2Right = -200:0.1:200;
yRange2Right = xRange2Right*mCoeff2Right + kCoeff2Right;
if ~isIdealCamera
    rad_fovHCamera = deg2rad(fovHCamera);
    mCoeff2Right_tmp = tan(pi/2-0.5*rad_fovHCamera)*sign(mCoeff2Right);
    if abs(mCoeff2Right_tmp) > abs(mCoeff2Right)
        mCoeff2Right = mCoeff2Right_tmp;
        yRange2Right = mCoeff2Right*(xRange2Right - cameraCenter(1)) + cameraCenter(2);
        kCoeff2Right = (xRange2Right(1)*yRange2Right(2) - xRange2Right(2)*yRange2Right(1))/(xRange2Right(1) - xRange2Right(2));
    end
end
plot(xRange2Right,yRange2Right,'k');

% Rette pasanti per gli estremi dell'elemento centrale
for j=1:floor(N)-1
    plot(XRange(j,:),YRange(j,:),'b')
end

% Calcolo degli estremi dell'elemento centrale
left_min_x = min(min(XRange));
left_min_y = min(YRange(XRange == left_min_x));
right_max_x = max(max(XRange));
right_max_y = min(YRange(find(XRange == right_max_x, 1 )));

% Retta passante per la camera e l'estremo sinistro dell'elemento centrale
y1 = left_min_y; x2 = cameraCenter(1);
x1 = left_min_x; y2 = cameraCenter(2);
mCoeff2Left_Central = (y1 - y2)/(x1 - x2);
kCoeff2Left_Central = (x1*y2 - x2*y1)/(x1 - x2);
xRange2Left_Central = -200:0.1:x1;
yRange2Left_Central = xRange2Left_Central*mCoeff2Left_Central + kCoeff2Left_Central;
plot(xRange2Left_Central,yRange2Left_Central,'k');

%Retta passante per la camera e l'estremo destro dell'elemento centrale
y1 = right_max_y; x2 = cameraCenter(1);
x1 = right_max_x; y2 = cameraCenter(2);
mCoeff2Right_Central = (y1 - y2)/(x1 - x2);
kCoeff2Right_Central = (x1*y2 - x2*y1)/(x1 - x2);
xRange2Right_Central = x1:0.1:200;
yRange2Right_Central = xRange2Right_Central.*mCoeff2Right_Central + kCoeff2Right_Central;
plot(xRange2Right_Central,yRange2Right_Central,'k');

% Retta passante per lo specchio sinistro
x1 = inPointLeft(1); x2 = exPointLeft(1);
y1 = -inPointLeft(2); y2 = -exPointLeft(2);
mCoeffLeft_Mirror = (y1 - y2)/(x1 - x2);
kCoeffLeft_Mirror = (x1*y2 - x2*y1)/(x1 - x2);

% Retta passante per lo specchio destro
y1 = -inPointRight(2); x2 = exPointRight(1);
x1 = inPointRight(1); y2 = -exPointRight(2);
mCoeffRight_Mirror = (y1 - y2)/(x1 - x2);
kCoeffRight_Mirror = (x1*y2 - x2*y1)/(x1 - x2);

vis_el_X = [];
vis_el_Y = [];
for i=1:floor(N)-1;
    for j=1:size(el_X,2)
        if (... % The point to be potentially visible must be located between the black lines (in previous plot) or below the central element, but "behind" the mirrors
            ( isOver(mCoeff2Left, kCoeff2Left, el_X(i,j), el_Y(i,j))&&isBelow(mCoeff2Left_Central, kCoeff2Left_Central, el_X(i,j), el_Y(i,j))||...
              isOver(mCoeff2Right, kCoeff2Right, el_X(i,j), el_Y(i,j))&&isBelow(mCoeff2Right_Central, kCoeff2Right_Central, el_X(i,j), el_Y(i,j)))||...
              isOver(mCoeffRight_Mirror, kCoeffRight_Mirror, el_X(i,j), el_Y(i,j))&&isOver(mCoeff2Right, kCoeff2Right, el_X(i,j), el_Y(i,j))&&el_Y(i,j)<right_max_y||...
              isOver(mCoeffLeft_Mirror, kCoeffLeft_Mirror, el_X(i,j), el_Y(i,j))&&isOver(mCoeff2Left, kCoeff2Left, el_X(i,j), el_Y(i,j))&& el_Y(i,j)<left_min_y)          

            vis_el_X(i,j) = el_X(i,j);
            vis_el_Y(i,j) = el_Y(i,j);
        end
    end
end

for i=1:size(vis_el_Y,1);
    tmp_x = vis_el_X(i,:);
    tmp_x(tmp_x==0) = [];
    tmp_y = vis_el_Y(i,:);
    tmp_y(tmp_y==0) = [];
    plot(tmp_x,tmp_y);
end

%% Drawing
num_p = 2; idx = []; num_c = 1; old_i = 0;
for i=1:size(v,2);
    [p1, p2] = getTangentLine(el_X(i,:),el_Y(i,:), cameraCenter);
    XX = [p1(1), p2(1)];
    YY = [p1(2), p2(2)];
    for k=1:2
        if (...  
        ( isOver(mCoeff2Left, kCoeff2Left, XX(k), YY(k))&&isBelow(mCoeff2Left_Central, kCoeff2Left_Central, XX(k), YY(k))||...
          isOver(mCoeff2Right, kCoeff2Right, XX(k), YY(k))&&isBelow(mCoeff2Right_Central, kCoeff2Right_Central, XX(k), YY(k)))||...
          isOver(mCoeffRight_Mirror, kCoeffRight_Mirror, XX(k), YY(k))&&isOver(mCoeff2Right, kCoeff2Right, XX(k), YY(k))&&YY(k)<right_max_y||...
          isOver(mCoeffLeft_Mirror, kCoeffLeft_Mirror, XX(k), YY(k))&&isOver(mCoeff2Left, kCoeff2Left, XX(k), YY(k))&& YY(k)<left_min_y)          
            if(~isBehind(el_X, el_Y, XX(k), YY(k), cameraCenter(1), cameraCenter(2), i))
                plot([XX(k) cameraCenter(1)], [YY(k) cameraCenter(2)]);
                num_p = num_p +1;
                if old_i ~= i 
                    idx = [idx,i];
                    num_c = num_c +1;
                    old_i = i;
                end
            end
        end
    end
end
scatter(nose(1),nose(2),'g');
%% Reflections counting
num_h = 1;
for i=1:size(idx,2);
    [XX, YY, rot] = rotateNose(nose, angleLeftMirror, angleRightMirror, mirrorsCenter, N, idx(i));
    if (...  
    ( isOver(mCoeff2Left, kCoeff2Left, XX, YY)&&isBelow(mCoeff2Left_Central, kCoeff2Left_Central, XX, YY)||...
      isOver(mCoeff2Right, kCoeff2Right, XX, YY)&&isBelow(mCoeff2Right_Central, kCoeff2Right_Central, XX, YY))||...
      isOver(mCoeffRight_Mirror, kCoeffRight_Mirror, XX, YY)&&isOver(mCoeff2Right, kCoeff2Right, XX, YY)&&YY<right_max_y||...
      isOver(mCoeffLeft_Mirror, kCoeffLeft_Mirror, XX, YY)&&isOver(mCoeff2Left, kCoeff2Left, XX, YY)&& YY<left_min_y)          
        if(~isBehind(el_X, el_Y, XX, YY, cameraCenter(1), cameraCenter(2), idx(i)))
            %if (~isNoseBehindHead(el_X(idx(i),:), el_Y(idx(i),:), XX, YY,  cameraCenter) && abs(rot) < 90)
            if (abs(rot) <= 89)
                num_h = num_h +1;
                scatter(XX,YY,'g');
            else
                scatter(XX,YY,'r');
            end
        end
    end
end
disp(['We can see always the frontal head but also ', num2str(num_h-1), ' Noses. There are ', num2str(num_c), ' heads in range. I can see ', num2str(num_p/2), ' partial heads.']);