%rappell fine tuning :



%###############    START      ######################
    %1 hyper para 
%Set initial default parameter
%default cut in spectrum dimension:
bottom_wl_treshold = 4;
top_wl_treshold = 30;
%default filter smoothing: medfilt2. todo: choose which default
%default nbr of clusters in weighted kmean:
nbr_of_cluster = 20;

%TODOOOOOO generlaize relevant cluster by tissue



%##############for one tissue:
tissue_name ='x40_r1_r2';
weighted_kmean_name = 'x24wkm_20_kmeans_and_D.mat';

crop_area_1 = [500 1000  (1000-500) (1500-1000)];
crop_area_2 = [1750 580  (2150-1750) (820-580)];
crop_area_3 = [ 1500 1000  (2000-1500) (1500-1000)];
crop_area_4 = [ 3400 1000  (3800-3400) (1400-1000)];
all = [];

seleted_area_list = [ crop_area_1 ; crop_area_2 ; crop_area_3 ; crop_area_4 ];
selected_area_name_list  = [ 'area_1';  'area_2';  'area_3'; 'area_4'];

selected_area = seleted_area_list(4,:);
%selected_area = all;


relevant_clusters = [0.21, 0.05, 0.37];  



%###########  select area we work on - if all the tissue m select area = all
%folder_path  = '/Users/lior/Desktop/Image & Analysis /cells_label/'+ string(tissue_name)+string('/')+selected_area_name+'/';


[data , lambda_mat , wl_mat ,rgb,  distance_map_list] = initialization_for_one_tissue(tissue_name, weighted_kmean_name, bottom_wl_treshold, top_wl_treshold, nbr_of_cluster, relevant_clusters); 


%############## get HnE
% if not working , le mettre plus bas.
[full_H, full_E] = get_HnE(wl_mat,lambda_mat);


%-----------------------------------------------------------------


%%%%% STEP AGADIR -  pr chaque area de tissue,  recupere les candidat_list
nbr_of_areas_by_tissue = 4;
%tissue_name_list = ['c161'; 'x40_r1_r2'];
!!!!!!!!!!!!!  tissue_name_list = ['x40_r1_r2'];
[nbr_of_tissues, ~] = size(tissue_name_list);
%tissue_area_aray_of_candidate_cells = [];
tissue_area_aray_of_candidate_cells = load('tissue_area_aray_of_candidate_cells.mat').tissue_area_aray_of_candidate_cells;

for t= 1:nbr_of_tissues
    %here insert all tissue calc including hnE
    for a=1: nbr_of_areas_by_tissue

        selected_area = seleted_area_list(a,:);
        selected_area_name =  selected_area_name_list(a,:);
        folder_path = pwd  + "/cells_label/" + string(tissue_name) + "/" +selected_area_name+ "/" ;

        [~ , s2_sa] = size(selected_area);
        if s2_sa == 0
            background = rgb;
            distance_map_list_a1 = distance_map_list(:,:,1);
            distance_map_list_a2 = distance_map_list(:,:,2);
            distance_map_list_a3 = distance_map_list(:,:,3);
            H = full_H;
            E = full_E;
        else
            %distance_map_list_full_a1 = distance_map_list(:,:,1);
            background = imcrop(rgb, selected_area);
            distance_map_list_a1 = imcrop(distance_map_list(:,:,1), selected_area);
            distance_map_list_a2 = imcrop(distance_map_list(:,:,2), selected_area);
            distance_map_list_a3 = imcrop(distance_map_list(:,:,3), selected_area);
            H = imcrop(full_H, selected_area);
            E = imcrop(full_E, selected_area);

        end

                %MSER
        ThresholdDelta =0.4;
        %,'MaxAreaVariation' ,.25);
        %'ThresholdDelta',3,
        
        
        %filter doublons parameters
                %MserRegions = regions_M_242_a2;
        %MserRegions = MergeM_a1;
        [bg_row,bg_col,~] = size(background);
        %max_aera = maxAereaRange;
        slide_by_raw = 200; 
        slide_by_col= 200;
        min_dist=20;
        plot_f_1 = false;
        bg=background;
        compress_percentage = .25;
        same_size_percentage = .85; %.6?
        min_range = 600;
        max_range = 1700;
        close_enough = 20; % a reduire en tant que premier filtre qui enleve les pures doublons. genre autour de 5-10


        [MserRegions_1,~] = detectMSERFeatures(distance_map_list_a1, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r1] =  filter_doublon_intra(MserRegions_1,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_2,~] = detectMSERFeatures(distance_map_list_a2, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r2] =  filter_doublon_intra(MserRegions_2,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_3,~] = detectMSERFeatures(distance_map_list_a3, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
 !!!!       [~, r3] =  filter_doublon_intra(MserRegions_3,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_H,~] = detectMSERFeatures(H, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, rH] =  filter_doublon_intra(MserRegions_H,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_E,~] = detectMSERFeatures(E, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, rE] =  filter_doublon_intra(MserRegions_E,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );


        concf_MSER_after_inter = MSERRegions([r1.PixelList;r2.PixelList;r3.PixelList; rH.PixelList; rE.PixelList]);
        

        %##############  Filter_1 : Doublons: removing closes/included ellipses

        [~,regions_after_inter_filter] = filter_doublon_inter(concf_MSER_after_inter ,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );
        %~~~~~~~~ ie F

        %figure('Name',selected_area_name);imshow(background); hold on;
        %hp = impixelinfo;
        %plot(regions_after_inter_filter,'showPixelList',false,'showEllipses',true);
        %hold off;
        
        
        tissue_area_aray_of_candidate_cells = [tissue_area_aray_of_candidate_cells ; {regions_after_inter_filter}, {full_H}, {full_E}, {folder_path}, {selected_area}  ] ;
    end
end

%save('tissue_area_aray_of_candidate_cells', 'tissue_area_aray_of_candidate_cells');
%tissue_area_aray_of_candidate_cells = load('tissue_area_aray_of_candidate_cells.mat').tissue_area_aray_of_candidate_cells;




%-----------------------------------------------------------------
%{

%%%%%%% STEP BARUCH : get confumat 

%[~,q] = size(tissue_area_aray_of_candidate_cells);  !!!!!!!!!!!! we dont implement it perfect. so q= t*a depends on t ! (a=4)



q = 1*4;
confumat = zeros(q, 21, 21,3);
for iou_treshold= 0:.05:1
    for HnE_quantile_treshold=0:.05:1
        for zone = 1:q
          %  if iou_treshold==0 && HnE_quantile_treshold==1 && zone==2
           %     dbg=1;
           % end
            candidates_regions = tissue_area_aray_of_candidate_cells{zone, 1};
            full_H = tissue_area_aray_of_candidate_cells{zone, 2};
            full_E = tissue_area_aray_of_candidate_cells{zone, 3};
            folder_p = tissue_area_aray_of_candidate_cells{zone, 4};   
            selected_area = tissue_area_aray_of_candidate_cells{zone, 5}; 

            gt_list = get_gt_list_from_folder_path (folder_p);

            regions_candidates_after_hne_filter_pixel_list = candidates_mser_regions_after_hne_filter(full_H, full_E,candidates_regions, HnE_quantile_treshold, selected_area);
            [s1,s2]= size(regions_candidates_after_hne_filter_pixel_list); 
            if s1 ~= 0
                 regions_candidates_after_hne_filter = MSERRegions(regions_candidates_after_hne_filter_pixel_list);      
            
                pred_list = get_pred_list_from_mser_regions(regions_candidates_after_hne_filter);


                [tp , fp, fn] = get_confusion_matrix_values(pred_list,gt_list , iou_treshold);  
            else
                tp = 0;
                fp = 0;
                [fn,~] = size(gt_list);
                
            end
            confumat(q,int16(iou_treshold *20)+1,int16(HnE_quantile_treshold *20)+1,1) = tp;
            confumat(q,int16(iou_treshold *20)+1,int16(HnE_quantile_treshold *20)+1,2) = fp;
            confumat(q,int16(iou_treshold *20)+1,int16(HnE_quantile_treshold *20)+1,3) = fn;
            
        end
    end
end
%-----------------------------------------------------------------
 

%%%%%% STEP CRYSTILLAZE
Precision_2d_mat = zeros(21,21);
Recall_2d_mat = zeros(21,21);

for iou_th_chiffre= 0:.05:1
    for HnE_quantile_th_grecque=0:.05:1
        tp_IoUchiffre_HnEgrq = sum(confumat( : ,int16(20 * iou_th_chiffre)+1 , int16(20*HnE_quantile_th_grecque)+1,1));
        fp_IoUchiffre_HnEgrq = sum(confumat( : ,int16(20 * iou_th_chiffre)+1 , int16(20*HnE_quantile_th_grecque)+1,2));
        fn_IoUchiffre_HnEgrq = sum(confumat( : ,int16(20 * iou_th_chiffre)+1 , int16(20*HnE_quantile_th_grecque)+1,3));

        if tp_IoUchiffre_HnEgrq==0 && fp_IoUchiffre_HnEgrq==0
            Precision_2d_mat(int16(20 * iou_th_chiffre)+1 ,  int16(20*HnE_quantile_th_grecque)+1)  = 1;
        else
            Precision_2d_mat(int16(20 * iou_th_chiffre)+1 ,  int16(20*HnE_quantile_th_grecque)+1)  = tp_IoUchiffre_HnEgrq / (tp_IoUchiffre_HnEgrq + fp_IoUchiffre_HnEgrq);
        end
        Recall_2d_mat(int16(20 * iou_th_chiffre)+1,  int16(20*HnE_quantile_th_grecque)+1) = tp_IoUchiffre_HnEgrq /( tp_IoUchiffre_HnEgrq + fn_IoUchiffre_HnEgrq);

    end
end

%-----------------------------------------------------------------


%%% STEP DORMS
precision_by_iou_th_list = [];
recall_by_iou_th_list = [];
Area_under_PR_curve_list_by_iou_th = [];

for iou_treshold= 0:.05:1
    precision_by_iou_th_list = [precision_by_iou_th_list; Precision_2d_mat(int16(20 * iou_treshold)+1 ,:)];
    recall_by_iou_th_list = [recall_by_iou_th_list; Recall_2d_mat(int16(20 * iou_treshold)+1,:)];
    Area_under_PR_curve = trapz(Precision_2d_mat(int16(20 * iou_treshold)+1,:),Recall_2d_mat(int16(20 * iou_treshold)+1,:));
    Area_under_PR_curve_list_by_iou_th = [Area_under_PR_curve_list_by_iou_th ;  Area_under_PR_curve];
end


%-----------------------------------------------------------------
            
%%%%%LOOP FOR VO+IZUALIZE I+DIFF CURVE AT ONCE 
%dici on peut deja fixer what's the best iou_treshold
lg =[]; 
for iou_treshold= 0:.05:1
lg= [ lg; 'IoU th =' + string(iou_treshold)];
end

plot(recall_by_iou_th_list(1,:) , precision_by_iou_th_list(1,:))
title('Precision-Recall graph for different IoU threshold')
hold on;
[ss,~] = size(precision_by_iou_th_list);
for i=2:ss
    plot(recall_by_iou_th_list(i,:) , precision_by_iou_th_list(i,:))
end
hold off; legend(lg, 'Location','southwest');

%-----------------------------------------------------------------
%}