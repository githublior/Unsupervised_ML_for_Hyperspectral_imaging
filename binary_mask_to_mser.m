function gt_cell_mser = binary_mask_to_mser(gt_cell_indix);
    [x,y] = size(gt_cell_indix);
    pixellist = [];
    for i=1:x
        for j=1:y
            if gt_cell_indix(i,j)
                pixellist = [pixellist; cast(j, 'int32' ) cast(i, 'int32')];
            end
        end
    end
    gt_cell_mser = MSERRegions({pixellist});
end