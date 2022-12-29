function [ gt_covering_rate, candidate_covering_rate] = compare_candidate_on_gt(candidate_cell, gt_cell);
    [nbr_pxl_in_candidate,~] = size(candidate_cell);
    match_counter =0;
    for i =1:nbr_pxl_in_candidate
        coord = candidate_cell(i,:);
        if gt_cell(coord(2)  , coord(1)) == true
            match_counter = match_counter +1;
        end        
    end
    nbr_pxl_in_gt =  sum(sum(gt_cell));
    gt_covering_rate = match_counter /nbr_pxl_in_gt;
    candidate_pixel_apport = match_counter/nbr_pxl_in_candidate;
    candidate_covering_rate = candidate_pixel_apport;
end