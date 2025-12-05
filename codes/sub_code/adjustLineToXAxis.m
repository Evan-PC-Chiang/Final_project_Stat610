function [x1_new, y1_new, x2_new, y2_new] = adjustLineToXAxis(x1, y1, x2, y2)
    % Initialize output with NaNs to indicate no valid segment if the line is skipped
    x1_new = NaN;
    y1_new = NaN;
    x2_new = NaN;
    y2_new = NaN;
    
    % Check if the line is totally above the x-axis
    if y1 > 0 && y2 > 0
        return;
    end
    
    % Check if the line crosses the y=0 axis
    if (y1 > 0 && y2 < 0) || (y1 < 0 && y2 > 0)
        % Find the intersection with the y = 0 line
        y_new = 0;
        t = (y_new - y1) / (y2 - y1);
        x_new = x1 + t * (x2 - x1);
        
        % Return the part of the line where y < 0
        if y1 < 0
            x1_new = x1;
            y1_new = y1;
            x2_new = x_new;
            y2_new = y_new;
        else
            x1_new = x_new;
            y1_new = y_new;
            x2_new = x2;
            y2_new = y2;
        end
    else
        % If the line is totally below the x-axis, return the original points
        x1_new = x1;
        y1_new = y1;
        x2_new = x2;
        y2_new = y2;
    end
end