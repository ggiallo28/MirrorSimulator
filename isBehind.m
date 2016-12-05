function bool_res = isBehind(all_head_X, all_head_Y, x, y, camera_x, camera_y, index)
    % y = (x*y1 - x*y2 + x1*y2 - x2*y1)/(x1 - x2);
    % y = x*(y1 - y2)/(x1 - x2) + (x1*y2 - x2*y1)/(x1 - x2);
    % m = (y1 - y2)/(x1 - x2)
    % k = (x1*y2 - x2*y1)/(x1 - x2)    
    x1 = camera_x; x2 = x;
    y1 = camera_y; y2 = y;
    mCoeff = (y1 - y2)/(x1 - x2);
    kCoeff = (x1*y2 - x2*y1)/(x1 - x2);
    DIFF_Y = [];
    for i=1:size(all_head_X,1)
        if(~(all_head_Y(i,:)>all_head_Y(index,:)))
            if(index ~= i)               
                for k=1:size(all_head_X,1)
                    if(~(all_head_Y(k,:)>all_head_Y(index,:)))
                        if(index ~= k)
                            yy = all_head_X(k,:).*mCoeff+kCoeff;
                            DIFF_Y = [DIFF_Y,abs(yy-all_head_Y(k,:))];
                        end
                    end
                end
            end
        end
    end
    idx = find(DIFF_Y<0.2);
    if(isempty(idx))
        bool_res = false;
    else
        bool_res = true;
    end


end