function dist = pointToRectangleDistance(x0, y0, x_min, y_min, x_max, y_max)
    % 计算点 (x0, y0) 到矩形的最小距离
    % 矩形的左下角 (x_min, y_min) 和右上角 (x_max, y_max)
    
    % 计算点到矩形的水平距离
    if x0 < x_min
        dx = x_min - x0;
    elseif x0 > x_max
        dx = x0 - x_max;
    else
        dx = 0;  % 点在矩形的水平边界内
    end
    
    % 计算点到矩形的垂直距离
    if y0 < y_min
        dy = y_min - y0;
    elseif y0 > y_max
        dy = y0 - y_max;
    else
        dy = 0;  % 点在矩形的垂直边界内
    end
    
    % 使用勾股定理计算点到矩形的最小距离
    dist = sqrt(dx^2 + dy^2);
end