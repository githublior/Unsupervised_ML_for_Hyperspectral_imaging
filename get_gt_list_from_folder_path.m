function gt_list = get_gt_list_from_folder_path (folder_path);
gt_path_list = get_cell_gt_list_from_folder_path(folder_path);
[s1,~] = size(gt_path_list);
gt_list = [];
for i=1:s1
    indiv_cell_gt = load(gt_path_list(i)).bw;
    gt_list = [gt_list; {indiv_cell_gt}];
end