function [tp , fp, fn] = get_confusion_matrix_values(pred_list, gt_list, IoU_Threshold);
tp =0;
fp =0;

[pred_s1,~] =  size(pred_list);
[gt_s1,~] =  size(gt_list);


%FLAGS                          % spy(sparse(current_gt_cell));
gt_flag = false([gt_s1 1]);
pred_flag = false([pred_s1 1]);


for pred_cell_index=1:pred_s1
    current_pred_cell = pred_list{pred_cell_index};
    for gt_cell_index = 1: gt_s1
        current_gt_cell = gt_list{gt_cell_index};

        compare = get_iou(current_pred_cell, current_gt_cell);
        if compare > IoU_Threshold
            %show_crossing(bg, [{current_pred_cell}], [{current_gt_cell}]);

            %current_pred_cell is TP. ie:
            tp = tp+1; 
            pred_flag(pred_cell_index) = true;
            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            %ATTENTION: si  gt_flag(gt_cell_index) etait deja positif,
            %alors doublon. ie 2 pred vise meme celll . que faire ds cette
            %situation ????
            gt_flag(gt_cell_index) = true;
            break
        end
        
    end
    
    if pred_flag(pred_cell_index)== false
        fp = fp +1;
    end
end
%les gt cells qui sont encore ac flag false:
fn = sum(~gt_flag);
end


%tp dimiue ac th augm         tp ie hit 
%fn augmente ac th augm       fn ie ceux que jai raté
%fp diminue ac conf_th augm       fp ie ceux que jai inventé