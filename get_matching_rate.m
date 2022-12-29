function [ matching_rate] = get_matching_rate(cell_pred, cell_gt);
    [l,~] = size(cell_pred);
    match =0;
    for i =1:l
        coord = cell_pred(i,:);
        if cell_gt(coord(2)  , coord(1)) == true
            match = match +1;
        end        
    end
    nbr_pxl_in_gt =  sum(sum(cell_gt));
    matching_rate = match /nbr_pxl_in_gt;
end