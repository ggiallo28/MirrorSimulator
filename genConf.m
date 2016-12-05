function Params = genConf()
    %% GENERATE PARAMETERS
    lengthMirrors = 65;         %cm
    mirror2Pivot = 2.5;         %cm Ci sta una struttura ad angolo retto che collega lo specchio all'ase di rotazione, abbiamo lungo l'asse dello specchio uno sfasamento di 2.5 cm
    offset2Mirror = 3;          %cm ed uno sfasamento di 3 cm in avanti.
    pivot2pivot = 28.5;         %cm Distanza tra i due perni
    angles = 30:2:80;
    cameradist = 120:10:220;
    Params = [];
    for count_i=1:size(angles,2)
        for count_k=1:size(cameradist,2)
            angleRightMirror = angles(count_i);
            angleLeftMirror = 180-angleRightMirror;
            lengthHypotenuse = sqrt(offset2Mirror^2 + mirror2Pivot^2); % Ipotenusa dello sfasamento
            angle2Pivot = rad2deg(asin(offset2Mirror/lengthHypotenuse));  %deg angolo in alto del triangolo sotto. Angolo fisso dovuto alla struttura (entità dell'offset)
            inPointLeft(1) =  lengthHypotenuse*cos(deg2rad(angleLeftMirror-angle2Pivot));                   %X 
            inPointLeft(2) =  lengthHypotenuse*sin(deg2rad(angleLeftMirror-angle2Pivot));                   %Y
            inPointRight(1) = pivot2pivot + lengthHypotenuse*cos(deg2rad(angleRightMirror+angle2Pivot));    %X, considero la traslazione lungo X per normalizzare i riferimenti
            inPointRight(2) = lengthHypotenuse*sin(deg2rad(angleRightMirror+angle2Pivot));                  %Y
            % Punti esterni specchi
            exPointLeft(1) = inPointLeft(1) + lengthMirrors*cos(deg2rad(angleLeftMirror));                  %X
            exPointLeft(2) = inPointLeft(2) + lengthMirrors*sin(deg2rad(angleLeftMirror));                  %Y
            exPointRight(1) = inPointRight(1) + lengthMirrors*cos(deg2rad(angleRightMirror));               %X,considero della traslazione lungo X
            exPointRight(2) = inPointRight(2) + lengthMirrors*sin(deg2rad(angleRightMirror));               %Y
            %% HEAD POSITION
            headPosx_range = 0.2:0.2:0.8;
            headPosy_range = 20:5:round(min([exPointRight(2),exPointRight(2)]));

            for count_j=1:size(headPosx_range,2)
                for count_l=1:size(headPosy_range,2)
                    Params = [Params; [angles(count_i), cameradist(count_k), headPosx_range(count_j)*pivot2pivot, headPosy_range(count_l)]];
                end
            end
        end
    end
end