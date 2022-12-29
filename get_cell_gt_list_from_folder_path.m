function gt_list = get_cell_gt_list_from_folder_path(folder_path);
    l =dir(folder_path);
    [q1,~]= size(l);
    gt_list = [];
    for lbl_idx=4:q1
        gt_list = [gt_list ; folder_path+l(lbl_idx).name];
    end

end