function [el_x, el_y, nose] = genHead(headPosx, headPosy, widthHead, lengthHead)
    el_x_1 = headPosx+0.5*widthHead*cos(-pi:0.01:pi);
    el_y_1 = headPosy+0.5*lengthHead*sin(-pi:0.01:pi);
    el_x_2 = headPosx+0.25*widthHead+0.1*widthHead*cos(-pi:0.01:pi);
    el_y_2 = headPosy+0.1*lengthHead*sin(-pi:0.01:pi);
    el_x = [el_x_1, el_x_2];
    el_y = [el_y_1, el_y_2];
    [nose_y, idx] = min(el_y);
    nose = [el_x(idx), nose_y];
    % TODO Raggiera per distinguere lato destro da sinistro
end

