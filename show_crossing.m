function show_crossing(bg, pred_list, gt_list)
    pred_cell_mser = MSERRegions(pred_list);
    
    
    figure('Name','gt & pred ');imshow(bg); hold on;
    hp = impixelinfo;
    plot(pred_cell_mser,'showPixelList',false,'showEllipses',true);

    [z,~] = size(gt_list);
    for c = 1:z
        gt_cell_index = gt_list{c,1};
        gt_cell_mser = binary_mask_to_mser(gt_cell_index);
        plot(gt_cell_mser,'showPixelList',true,'showEllipses',false);
    end
end