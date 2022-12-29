function [fpr, tpr] = cell_segm_get_score_(struct_path,pred_all, gt_all_cell_union);
% calculate fpr tpr pf the prediction, wi default parameter like percentage
% of matching between pred and gt, and...?

%tpr = tp/p =  tp/(tp+fn)
%fpr = fp /n = fp/ (fp+tn)

% check cbn (qtite) pred_i pxls are in gt_mat==1
%for pixel=all elmt of pred_i:
%   if gt_mat(pixel)==1:
        %match +=1
    %else: mismatch +=1
%end




%end
%didnt_catched = len(gt_mat)- match
[a,b]  = size(pic);
 all_pred_as_false = a*b - len(pred)
 all_gt_as_false = a*b - sum(gt_all_cell_union)
tp = match
fn = didnt_catched  !FAUXXXXX
fp = mismatch
tn = len( ???)

nmb_pxl)in_img = a*b
??? =  (a*b - len(pred))- 