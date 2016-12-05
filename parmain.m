clc; clear all; close all;
it = 0;
set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"
warning off
%% Costants
lengthMirrors = 65;         %cm
mirror2Pivot = 2.5;         %cm Ci sta una struttura ad angolo retto che collega lo specchio all'ase di rotazione, abbiamo lungo l'asse dello specchio uno sfasamento di 2.5 cm
offset2Mirror = 3;          %cm ed uno sfasamento di 3 cm in avanti.
pivot2pivot = 28.5;         %cm Distanza tra i due perni
lengthHead = 15;            %cm % vertical 2*radius
widthHead = 10;             %cm % horizontal 2*radius
%% Parameters
angles = 30:80;
cameradist = 120:5:220;
Params = genConf();
Results = cell(size(Params,1),1);
parfor count = 1:size(Params,1)
    if exist(['D:/Plot/plot',num2str(count),'.png'], 'file')
        toplot = false;
    else
        toplot = true;
    end
    toplot1 = false;
                
    angleRightMirror = Params(count,1);
    angleLeftMirror = 180-angleRightMirror;      %deg
    % distanceCamera = 220;       %cm
    distanceCamera = Params(count,2);       %cm
    fovHCamera = 50;            %deg
    isIdealCamera = true;
    %% Logic
    lengthPivotMirrors = lengthMirrors+mirror2Pivot; % Lunghezza dello specchio più offset dovuto al fatto che lo specchio non finisce dove è posto l'asse di rotazione (pivot).
    angleBetweenMirror = abs(angleRightMirror-angleLeftMirror); % Angolo tra gli specchi
    lengthHypotenuse = sqrt(offset2Mirror^2 + mirror2Pivot^2); % Ipotenusa dello sfasamento
    angle2Pivot = rad2deg(asin(offset2Mirror/lengthHypotenuse));  %deg angolo in alto del triangolo sotto. Angolo fisso dovuto alla struttura (entità dell'offset)
    N = 360/angleBetweenMirror; % FORMULA 1: Dice quanti riflessi vengono prodotti, quelli visibili saranno un sottoinsieme. Da usare come verifica.
    % Proietto l'ipotenusa che si forma sull'asse X e Y in modo da avere le coordinate del punto in basso a destra. Il sistema
    % di riferimento è posto sul perno sinisto.
    %   ---------------->                          -------------------|->
    %   |\                                                           /|
    %   | \                                                         / |   %%  Standardizzazione dei riferimenti: centro sul perno sinistro
    %   |__\ <- Punto                                    Punto ->  /__|
    %       |                                                     |
    % Calcolo posizione punti interni specchi (estremi inferiori)
    inPointLeft = []; inPointLeft(1) =  lengthHypotenuse*cos(deg2rad(angleLeftMirror-angle2Pivot));                   %X 
    inPointLeft(2) =  lengthHypotenuse*sin(deg2rad(angleLeftMirror-angle2Pivot));                   %Y
    inPointRight = []; inPointRight(1) = pivot2pivot + lengthHypotenuse*cos(deg2rad(angleRightMirror+angle2Pivot));    %X, considero la traslazione lungo X per normalizzare i riferimenti
    inPointRight(2) = lengthHypotenuse*sin(deg2rad(angleRightMirror+angle2Pivot));                  %Y
    distanceBetweenMirror = []; distanceBetweenMirror(1) = pdist([inPointLeft;inPointRight],'euclidean'); %shortest
    % Punti esterni specchi
    exPointLeft = []; exPointLeft(1) = inPointLeft(1) + lengthMirrors*cos(deg2rad(angleLeftMirror));                  %X
    exPointLeft(2) = inPointLeft(2) + lengthMirrors*sin(deg2rad(angleLeftMirror));                  %Y
    exPointRight = []; exPointRight(1) = inPointRight(1) + lengthMirrors*cos(deg2rad(angleRightMirror));               %X,considero della traslazione lungo X
    exPointRight(2) = inPointRight(2) + lengthMirrors*sin(deg2rad(angleRightMirror));               %Y
   
    headPosx = Params(count,3); %cm % x0,y0 ellipse centre coordinates
    headPosy = Params(count,4);             %cm
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
    mirrorsCenter = []; [mirrorsCenter(1),mirrorsCenter(2)] = inters2rette(mLeft,kLeft,mRight,kRight);
    % Sto ipotizzando la camera al centro (X) e posta ad una distanza (Y) definita precedentemente
    cameraCenter = []; cameraCenter(1) = 0.5*pivot2pivot; cameraCenter(2) = -distanceCamera; 
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
    if toplot1
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
        if toplot1
            plot(XRange(j,:),YRange(j,:))
        end
        new_center_x = mean(XRange(j,:));
        new_center_y = mean(YRange(j,:));
        distROT = sqrt((mirrorsCenter(1)-mean(XRange(j,:)))^2 + (mirrorsCenter(2)-mean(YRange(j,:)))^2);
        assert(abs(distOR-distROT)<exp(-10));
    end
    if toplot
        hold off;
        %% Plot delle teste potenzialmente visibili
        fig = figure; hold on,axis([-max(cameradist) max(cameradist) -max(cameradist) max(cameradist)])
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
    end

    % Retta passante per il centro della camera e l'estremo esterno dello specchio sinistro
    y1 = -exPointLeft(2); x2 = cameraCenter(1);
    x1 = exPointLeft(1); y2 = cameraCenter(2);
    mCoeff2Left = (y1 - y2)/(x1 - x2);
    kCoeff2Left = (x1*y2 - x2*y1)/(x1 - x2);
    xRange2Left = -200:0.1:200;
    yRange2Left = xRange2Left*mCoeff2Left + kCoeff2Left;
    if ~isIdealCamera
        rad_fovHCamera = deg2rad(fovHCamera);
        mCoeff2Left = tan(pi/2-0.5*rad_fovHCamera)*sign(mCoeff2Left);
        yRange2Left = mCoeff2Left*(xRange2Left - cameraCenter(1)) + cameraCenter(2);  
        kCoeff2Left = (xRange2Left(1)*yRange2Left(2) - xRange2Left(2)*yRange2Left(1))/(xRange2Left(1) - xRange2Left(2));
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
        mCoeff2Right = tan(pi/2-0.5*rad_fovHCamera)*sign(mCoeff2Right);
        yRange2Right = mCoeff2Right*(xRange2Right - cameraCenter(1)) + cameraCenter(2);
        kCoeff2Right = (xRange2Right(1)*yRange2Right(2) - xRange2Right(2)*yRange2Right(1))/(xRange2Right(1) - xRange2Right(2));
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
    num_p = 2;
    num_c = 1;
    old_i = 0;
    idx = [];
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
                    if toplot
                        plot([XX(k) cameraCenter(1)], [YY(k) cameraCenter(2)]);
                    end
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
                if (~isNoseBehindHead(el_X(idx(i),:), el_Y(idx(i),:), XX, YY,  cameraCenter) && abs(rot) < 90)
                    num_h = num_h +1
                    scatter(XX,YY,'g');
                else
                    scatter(XX,YY,'r');
                end
            end
        end
    end
    % Params = [Params; [angles(count_i), cameradist(count_k), headPosx_range(count_j)*pivot2pivot, headPosy_range(count_l)]];
    Results(count) = {[Params(count,:),num_c,num_p/2,num_h]};
    if toplot
        print(fig,['D:/Plot/plot',num2str(count)],'-dpng');
        close all;
    end
    disp(['Iteration: ',num2str(count),' N° Com: ', num2str(num_c), ' N° Part: ',num2str(num_p/2), ' Angle: ' , num2str(Params(count,1))]);
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

