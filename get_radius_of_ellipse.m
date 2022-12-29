function rd = get_radius_of_ellipse(roi);
    list = roi.Vertices;
    [a1,a2] = size(list)

    low_y=  9999999999999;
    high_y = -1;
    for runner=1:a1
        low_y = min(low_y,list(runner, 1));
        high_y = max(high_y, list(runner,1));
    end
    
    
    low_x=  9999999999999;
    high_x = -1; 
    for runner=1:a2
        low_x = min(low_x,list(runner, 1));
        high_x = max(high_x, list(runner,1));
    end
    
    y_rd = (high_y - low_y)/2;
    x_rd = (high_x - low_x)/2;
    
    rd = max(y_rd, x_rd);

end