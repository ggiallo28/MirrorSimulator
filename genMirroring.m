function [v] = genMirroring( angleLeftMirror, angleRightMirror, anglePoint, N)
    h2m(1) = abs(angleLeftMirror-anglePoint); % Angle 2 Left -1
    h2m(2) = abs(angleRightMirror-anglePoint); % Angle 2 Right 1
    %% Compute Rotations, Definisco sfasamenti in angoli tra le teste
    % v(1,:) = angoli di rotazione
    % v(2,:) = punto da ruotare
    % v(3,:) = direzione del flip
    isLeft = true; l = 2; r = 1; v = zeros(2,floor(N)-1);
    v(1,1) = -2*h2m(1); v(1,2) = 2*h2m(2);
    for i=3:floor(N)-1
        if(isLeft)
            v(1,i) = v(1,i-2) - 2*h2m(l);
            if l==1
                l = 2;
            else
                l = 1;
            end
        else
            v(1,i) = v(1,i-2) + 2*h2m(r);
            if r==1
                r = 2;
            else
                r = 1;
            end
        end
        isLeft = ~isLeft;
    end
    
    k = 1;
    for i=1:floor(N)-1  
       v(2,i) = k;
       if mod(i,2) == 0 && i>0
           k = ~k;
       end
    end    
end


