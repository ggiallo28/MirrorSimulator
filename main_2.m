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
angleRightMirror = 80;                       %deg
angleLeftMirror = 180-angleRightMirror;      %deg
distanceCamera = 100;                        %cm
fovHCamera = 50;                             %deg
isIdealCamera = true;                        %the value fovHCamera is used only if the value isIdealCamera is true, otherwise we consider a camera with variable field of view
%% HEAD MODEL
lengthHead = 15;            %cm % vertical 2*radius
widthHead = 10;             %cm % horizontal 2*radius
headPosx = 0.5*pivot2pivot;                  %cm
headPosy = 12;                               %cm
%% Logic
angleBetweenMirror = abs(angleRightMirror-angleLeftMirror);     % Angle betwwen the two mirrors
N = 360/angleBetweenMirror; % FORMULA 1: Dice quanti riflessi vengono prodotti, quelli visibili saranno un sottoinsieme. Da usare come verifica.
% Proietto l'ipotenusa che si forma sull'asse X e Y in modo da avere le coordinate del punto in basso a destra. Il sistema
% di riferimento è posto sul perno sinisto.
%   ---------------->                          -------------------|->
%   |\                                                           /|
%   | \                                                         / |   %%  Standardizzazione dei riferimenti: centro sul perno sinistro
%   |__\ <- Punto                                    Punto ->  /__|
%       |                                                     |
% Calcolo posizione punti interni specchi (estremi inferiori)
inPointLeft(1) =  lengthHypotenuse*cos(deg2rad(angleLeftMirror-angle2Pivot));                   %X 
inPointLeft(2) =  lengthHypotenuse*sin(deg2rad(angleLeftMirror-angle2Pivot));                   %Y
inPointRight(1) = pivot2pivot + lengthHypotenuse*cos(deg2rad(angleRightMirror+angle2Pivot));    %X, considero la traslazione lungo X per normalizzare i riferimenti
inPointRight(2) = lengthHypotenuse*sin(deg2rad(angleRightMirror+angle2Pivot));                  %Y
distanceBetweenMirror(1) = pdist([inPointLeft;inPointRight],'euclidean'); %shortest
% Punti esterni specchi
exPointLeft(1) = inPointLeft(1) + lengthMirrors*cos(deg2rad(angleLeftMirror));                  %X
exPointLeft(2) = inPointLeft(2) + lengthMirrors*sin(deg2rad(angleLeftMirror));                  %Y
exPointRight(1) = inPointRight(1) + lengthMirrors*cos(deg2rad(angleRightMirror));               %X,considero della traslazione lungo X
exPointRight(2) = inPointRight(2) + lengthMirrors*sin(deg2rad(angleRightMirror));               %Y
distanceBetweenMirror(2) = pdist([exPointLeft;exPointRight],'euclidean');
% Coefficienti angolari e costante delle rette passanti per i due specchi, mi serve per calcolare la posizione del centro di rotazione
% y = (x*y1 - x*y2 + x1*y2 - x2*y1)/(x1 - x2);
% y = x*(y1 - y2)/(x1 - x2) + (x1*y2 - x2*y1)/(x1 - x2);
% m = (y1 - y2)/(x1 - x2)
% k = (x1*y2 - x2*y1)/(x1 - x2)
mLeft = ((-inPointLeft(2)) - (-exPointLeft(2)))/(inPointLeft(1) - exPointLeft(1));
kLeft = (inPointLeft(1)*(-exPointLeft(2)) - exPointLeft(1)*(-inPointLeft(2)))/(inPointLeft(1) - exPointLeft(1));
mRight = ((-inPointRight(2)) - (-exPointRight(2)))/(inPointRight(1) - exPointRight(1));
kRight = (inPointRight(1)*(-exPointRight(2)) - exPointRight(1)*(-inPointRight(2)))/(inPointRight(1) - exPointRight(1));
[mirrorsCenter(1),mirrorsCenter(2)] = inters2rette(mLeft,kLeft,mRight,kRight);
% %% HEAD POSITION
% xH = -200:200;
% m_dist = 64;
% xoff = sqrt(m_dist/(1-cos(deg2rad(0.5*angleBetweenMirror))^2));
% yHLeft = mLeft*(xH - mirrorsCenter(1)) + mirrorsCenter(2) - xoff;  
% kHLeft = (xH(1)*yHLeft(2) - xH(2)*yHLeft(1))/(xH(1) - xH(2));
% 
% yHRight = mRight*(xH - mirrorsCenter(1)) + mirrorsCenter(2) - xoff;  
% kHRight = (xH(1)*yHRight(2) - xH(2)*yHRight(1))/(xH(1) - xH(2));
% 
% aa = xoff*cos(deg2rad(0.5*angleBetweenMirror));
% bb = sqrt(xoff*xoff-aa*aa)
% %headPosy = 16; headPosx = (headPosy-kHRight)/mRight; headPosx =
% %(headPosy-kHLeft)/mLeft;
% % plot(xH,yHLeft) plot(xH,yHRight)

% Sto ipotizzando la camera al centro (X) e posta ad una distanza (Y) definita precedentemente
cameraCenter(1) = 0.5*pivot2pivot; cameraCenter(2) = -distanceCamera;
if ~isIdealCamera
    cameraCenter(2) = -(exPointRight(2) +0.5*distanceBetweenMirror(2)*sqrt(1/cos(deg2rad(180-90-0.5*fovHCamera))^2-1));
    distanceCamera = abs(cameraCenter(2));
end
% Angolo della retta passante per il centro degli specchi e il centro della testa, serve per capire di quanto ruotare
% Deommentare per allineare la testa e la camera con il centro di rotazione degli specchi
% headPosx = mirrorsCenter(1);
% cameraCenter(1) = mirrorsCenter(1);
angleHead = rad2deg(atan2((-headPosy-mirrorsCenter(2)),(headPosx-mirrorsCenter(1))));
if angleHead<0
    angleHead = -angleHead;
end
vH = genMirroring(angleLeftMirror, angleRightMirror, angleHead, N);

% [el_x, el_y] = genHead(headPosx, -headPosy, widthHead, lengthHead);
% v = genMirroring(angleLeftMirror, angleRightMirror, angleHead, N);
% el_X = zeros(size(v,2),size(el_x,2));
% el_Y = zeros(size(v,2),size(el_y,2));
% distOR = sqrt((mirrorsCenter(1)-mean(el_x))^2 + (mirrorsCenter(2)-mean(el_y))^2);
% %% Flip per effetto mirroring
% for i=1:size(v,2);
%     el_x_flip = el_x;
%     el_y_flip = el_y;
%     if v(2,i) ~= 0
%         [el_x_flip,el_y_flip] = rotate(el_x_flip,el_y_flip,headPosx,-headPosy,180,'y');
%     end
%     [el_X(i,:),el_Y(i,:)] = rotate(el_x_flip,el_y_flip,mirrorsCenter(1),mirrorsCenter(2),v(1,i),'z');
%     new_center_x = mean(el_X(i,:));
%     new_center_y = mean(el_Y(i,:));
%     distROT = sqrt((mirrorsCenter(1)-mean(el_X(i,:)))^2 + (mirrorsCenter(2)-mean(el_Y(i,:)))^2);
%     assert(abs(distOR-distROT)<exp(-10));
% end

[el_x, el_y, nose] = genHead(headPosx, -headPosy, widthHead, lengthHead);
el_X = zeros(floor(N)-1,size(el_x,2));
el_Y = zeros(floor(N)-1,size(el_y,2));
distOR = sqrt((mirrorsCenter(1)-mean(el_x))^2 + (mirrorsCenter(2)-mean(el_y))^2);
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
for j=1:floor(N)-1
    new_center_x = mean(el_X(j,:));
    new_center_y = mean(el_Y(j,:));
    distROT = sqrt((mirrorsCenter(1)-mean(el_X(j,:)))^2 + (mirrorsCenter(2)-mean(el_Y(j,:)))^2);
    assert(abs(distOR-distROT)<exp(-10));
end

%% Plot Inverti le Y. plot([x1 x2], [y1 y2])
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

%% Disegnare elemento centrale
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
%% Plot delle teste potenzialmente visibili
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

% Retta passante per il centro della camera e l'estremo esterno dello specchio sinistro
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

% Retta passante per il centro della camera e l'estremo esterno dello specchio destro
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
        if (... % Il punto per essere potenzialmente visibile deve trovarsi tra le linee nere oppure al di sotto dell'elemento centrale ma "dietro" gli specchi
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

%% Plot delle proiezioni realmente visibili
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
%% Conta numero riflessi
num_h = 1;
for i=1:size(idx,2);
    [XX, YY, rot] = rotateNose(nose, angleLeftMirror, angleRightMirror, mirrorsCenter, N, idx(i));
    if (...  
    ( isOver(mCoeff2Left, kCoeff2Left, XX, YY)&&isBelow(mCoeff2Left_Central, kCoeff2Left_Central, XX, YY)||...
      isOver(mCoeff2Right, kCoeff2Right, XX, YY)&&isBelow(mCoeff2Right_Central, kCoeff2Right_Central, XX, YY))||...
      isOver(mCoeffRight_Mirror, kCoeffRight_Mirror, XX, YY)&&isOver(mCoeff2Right, kCoeff2Right, XX, YY)&&YY<right_max_y||...
      isOver(mCoeffLeft_Mirror, kCoeffLeft_Mirror, XX, YY)&&isOver(mCoeff2Left, kCoeff2Left, XX, YY)&& YY<left_min_y)          
        if(~isBehind(el_X, el_Y, XX, YY, cameraCenter(1), cameraCenter(2), idx(i)))
            rot
            if (~isNoseBehindHead(el_X(idx(i),:), el_Y(idx(i),:), XX, YY,  cameraCenter) && abs(rot) < 90)
                num_h = num_h +1
                scatter(XX,YY,'g');
            else
                scatter(XX,YY,'r');
            end
        end
    end
end
disp(['We can see always the frontal head but also ', num2str(num_h-1), ' Nose. There are ', num2str(num_c), ' heads in range. I can see ', num2str(num_p/2), ' partial head.']);
%% Evaluation
xX = -200:200;
m_dist = 64;
xoff = sqrt(m_dist/(1-cos(deg2rad(0.5*angleBetweenMirror))^2));
yHLeft = mLeft*(xX - mirrorsCenter(1)) + mirrorsCenter(2) - xoff;  
kHLeft = (xX(1)*yHLeft(2) - xX(2)*yHLeft(1))/(xX(1) - xX(2));

yHRight = mRight*(xX - mirrorsCenter(1)) + mirrorsCenter(2) - xoff;  
kHRight = (xX(1)*yHRight(2) - xX(2)*yHRight(1))/(xX(1) - xX(2));

aa = xoff*cos(deg2rad(0.5*angleBetweenMirror));
bb = sqrt(xoff*xoff-aa*aa);
if isIdealCamera
    fovHCamera = rad2deg(atan(abs((mCoeff2Left-mCoeff2Right)/(1+mCoeff2Left*mCoeff2Right))));
    [diamondAngle, idx] = min([abs(fovHCamera + 2*angleBetweenMirror),360-(fovHCamera + 2*angleBetweenMirror)]);
else
    [diamondAngle, idx] = min([fovHCamera + 2*angleBetweenMirror,360-(fovHCamera + 2*angleBetweenMirror)]);
end
delta = deg2rad(0.5*(180-diamondAngle));
%if nose(2) > min(-exPointRight(2), -exPointLeft(2))
if idx == 1
    mLL = tan(delta);
    mRR = tan(-delta);
else
    mLL = tan(-delta);
    mRR = tan(+delta);
end
    yLL = mLL*(xX-exPointLeft(1)) - exPointLeft(2);
    yRR = mRR*(xX-exPointRight(1)) - exPointRight(2);
    kLL = (xX(1)*yLL(2) - xX(2)*yLL(1))/(xX(1) - xX(2));
    kRR = (xX(1)*yRR(2) - xX(2)*yRR(1))/(xX(1) - xX(2));
    [xLim, yLim] = inters2rette(mLL,kLL,mRR,kRR);
if idx == 1
    xRR = exPointLeft(1):xLim;
    xLL = xLim:exPointRight(1);
else
    xLL = exPointLeft(1):xLim;
    xRR = xLim:exPointRight(1);
end   
    yLL = mLL*(xLL-exPointLeft(1)) - exPointLeft(2);
    yRR = mRR*(xRR-exPointRight(1)) - exPointRight(2);
    plot(xRR,yRR)
    plot(xLL,yLL)

    if angleBetweenMirror<45
        Limits = genMirroring(angleLeftMirror, angleRightMirror, angleHead, N);
        mRadiusLeft = tan(deg2rad(-90+Limits(1,3)));
        mRadiusRight = tan(deg2rad(-90+Limits(1,4)));
        yRadiusL = mRadiusLeft*(xX-mirrorsCenter(1)) + mirrorsCenter(2);
        yRadiusR = mRadiusRight*(xX-mirrorsCenter(1)) + mirrorsCenter(2);
        kRadiusLeft = (xX(1)*yRadiusL(2) - xX(2)*yRadiusL(1))/(xX(1) - xX(2));
        kRadiusRight = (xX(1)*yRadiusR(2) - xX(2)*yRadiusR(1))/(xX(1) - xX(2));
        [xLimRadiusL, yLimRadiusL] = inters2rette(mRadiusLeft,kRadiusLeft,mCoeff2Left,kCoeff2Left);
        [xLimRadiusR, yLimRadiusR] = inters2rette(mRadiusRight,kRadiusRight,mCoeff2Right,kCoeff2Right);
        Radius(1) = pdist([xLimRadiusL,yLimRadiusL;mirrorsCenter],'euclidean');
        Radius(2) = pdist([xLimRadiusR,yLimRadiusR;mirrorsCenter],'euclidean');
        mRadiusLeft = tan(atan(mCoeff2Left) - pi/2);
        mRadiusRight = tan(atan(mCoeff2Right) + pi/2);
        yRadiusL = mRadiusLeft*(xX-mirrorsCenter(1)) + mirrorsCenter(2);
        yRadiusR = mRadiusRight*(xX-mirrorsCenter(1)) + mirrorsCenter(2);
        kRadiusLeft = (xX(1)*yRadiusL(2) - xX(2)*yRadiusL(1))/(xX(1) - xX(2));
        kRadiusRight = (xX(1)*yRadiusR(2) - xX(2)*yRadiusR(1))/(xX(1) - xX(2));
        [xLimRadiusL, yLimRadiusL] = inters2rette(mRadiusLeft,kRadiusLeft,mCoeff2Left,kCoeff2Left);
        [xLimRadiusR, yLimRadiusR] = inters2rette(mRadiusRight,kRadiusRight,mCoeff2Right,kCoeff2Right);
        Radius(3) = pdist([xLimRadiusL,yLimRadiusL;mirrorsCenter],'euclidean');
        Radius(4) = pdist([xLimRadiusR,yLimRadiusR;mirrorsCenter],'euclidean');
        yLim5 = mirrorsCenter(2)-min(Radius)-0.5*lengthHead;

        yLL5 = yLL-yLim+yLim5;
        yRR5 = yRR-yLim+yLim5;
    %     plot(xRR,yRR5)
    %     plot(xLL,yLL5)
        xCircle = mirrorsCenter(1)+min(Radius)*cos(-pi:0.01:pi);
        yCircle = mirrorsCenter(2)+min(Radius)*sin(-pi:0.01:pi);

        plot(xCircle,yCircle)

        mRadiusLeft = tan(deg2rad(-90-2*(45-angleBetweenMirror)));
        mRadiusRight = tan(deg2rad(-90+2*(45-angleBetweenMirror)));
        yRadiusL = mRadiusLeft*(xX-mirrorsCenter(1)) + mirrorsCenter(2);
        yRadiusR = mRadiusRight*(xX-mirrorsCenter(1)) + mirrorsCenter(2);
        plot(xX,yRadiusL);
         plot(xX,yRadiusR);
    end
%% Test se abbiamo intersezione
% radius = pdist([mirrorsCenter;[headPosx, -headPosy]],'euclidean');
% %circle(mirrorsCenter(1),mirrorsCenter(2),radius);
% % Se abbiamo intersezione basta trovare l'intersezione tra le rette che passano per gli estremi degli specchi e il centro camera, in questo modo otteniamo il settore 
% % angolare dove le teste sono visibili.
% P1 = [-35.415,13.75]; % Punto a sinistra
% P2 = [63.73,13.76]; % Punto a destra
% PC = [14.1200,-20]; % Punto centrale
% d1 = pdist([P1;PC],'euclidean');
% d2 = pdist([P2;PC],'euclidean');
% a1 = asin(d1/(2*radius)); % settore angolare sinistro, da considerarsi positivo data la convenzione vH e il sistema di riferimento usato
% a2 = asin(d2/(2*radius)); % settore angolare destro, da considerarsi negativo data la convenzione vH e il sistema di riferimento usato
% % Il sistema di riferimento ha X -> Y verso il basso e Z entrante: è centrato dove si trova il pivot sinistro.
% 
% % Se non abbiamo intersezione (delta negativo) bisogna considerare il punto tangente
% % EQ1 = y - yd - m*(x - xd)
% % EQ2 = -r + (x-xc)^2 + (y-yc)^2 % Sostitusico la nella 2 con la 1
% syms r x xc yc xd yd r m 
% EQ = -r + (x-xc)^2 + (yd + m*(x - xd)-yc)^2;
% c = -r + (0-xc)^2 + (yd + m*(0 - xd)-yc)^2;
% a = m^2 + 1;
% b = subs(EQ - c - x*x*a,x,1);
% solve(b*b-4*a*c,m)
% Il punto tangente è caratterizzato dall'avere un delta nullo.
% trovati i settori angolari saranno visibili le sole teste che avranno angolo compreso tra il minimo e il massimo settore angolare
%% TODO calcolare indice di copertura

