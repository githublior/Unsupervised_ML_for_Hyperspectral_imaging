function [iou] = get_iou(cell_pred, cell_gt);
    [l,b] = size(cell_pred);
    inters =0;
    union_mat = false([l b]);

    for i =1:l
        coord = cell_pred(i,:);
        if cell_gt(coord(2)  , coord(1)) == true
            inters = inters +1;
        end
        union_mat(coord(2)  , coord(1)) = true;
    end
    uni = sum(sum(union_mat));
    
    iou= inters/uni;
end