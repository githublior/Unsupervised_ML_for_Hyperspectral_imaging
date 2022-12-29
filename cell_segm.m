
%[data, hne ] = preprocess_tissue(tissue_struct)


%{
tissue_name = "3.rawstack";
imaging = 'x20';
data_path = fullfile(pwd, tissue_name + ".mat");

data = load(data_path);
rgb = data.RGB;
wl_mat = data.spec;
lambda_mat = data.lambda;
[full_H, full_E] = get_HnE(wl_mat,lambda_mat);
hne.H = full_H;
hne.E = full_E;
dbg=1;
[regions_candidates_after_hne_filter_pixel_list] = cell_segm(data , hne, imaging, dbg)
%}


function [pred_list] = cell_segm(data , hne, imaging, dbg)
%data: struct containing data.spec : matrix of the wavelength
%spectrum,  data.lambda : the frequency of each wavelength channel,
%data.rgb : rgb visualization of the tissue. the data is supposed to be
%already preprocessed (preprocess_tissue.m)
%tissue_H_map: Hematoxylin matrix of the tissue
%tissue_E_map: Eosin matrix of the tissue




rgb = data.RGB;
wl_mat = data.spec;
lambda_mat = data.lambda;

H = hne.H;
E = hne.E;




%%%%%%%   Initial_Parameter   %%%%%%%%


%default filter smoothing: medfilt2. todo: choose which default
%default nbr of clusters in weighted kmean:
warning('off');
nbr_of_cluster = 20;
plot_f_1 = false;
HnE_quantile_treshold = .75;
compress_percentage = .25;
same_size_percentage = .85;
slide_by_raw = 200; 
slide_by_col= 200;
close_enough = 20;
min_dist=20;
[bg_row,bg_col,~] = size(rgb);
ThresholdDelta =0.4;







%weighted KMeans
%or load weighted kmeans 
if dbg==1
    wkm= load('3.rawstack_smoothed_and_wl_shaped_20_kmeans_and_D.mat'); 
    W_kmean_mat_output = wkm.kmean_mat_output;
    D = wkm.D;
else

    % Weighted KMeans computation
    fprintf('computing Weihted KMeans \n');
    input_spectral_img = wl_mat; %wl_mat
    %input_spectral_img = imcrop(input_spectral_img,[2761 686  (3781-2761) (1409-686)]);
    number_of_clusters = nbr_of_cluster;
    save_kmean_img=0; %1 if want to save img
    name_of_input_img = 'specific_tissue_wtvr';
    [W_kmean_mat_output, D]=weighted_kmeans_algo(input_spectral_img,number_of_clusters, save_kmean_img, name_of_input_img,pwd);
    fprintf('Weihted KMeans terminated\n');
end



%cell size selection 

if imaging == 'x40'
    mser_min_range = 200;
    mser_max_range = 8000;
else
    mser_min_range = 100;
    mser_max_range = 4000; 
end
%mser_max_range = get_cell_size_from_user(rgb, 'large');
%mser_min_range = get_cell_size_from_user(rgb, 'small');


% User's cluster selection
relevant_cluster = get_user_cluster_selection(W_kmean_mat_output);


distance_map_list = postprocess_on_cluster_selection(relevant_cluster, nbr_of_cluster,  D, W_kmean_mat_output);


        
[~,~, dmp_size]  = size(distance_map_list);
lmm = [];
for i=1:dmp_size

        fprintf('computing for cluster %d \n', i);
        [MserRegions_i,~] = detectMSERFeatures(distance_map_list(:,:,i), 'RegionAreaRange',[mser_min_range ,mser_max_range],'ThresholdDelta',ThresholdDelta);
        [~, ri] =  filter_doublon_intra(MserRegions_i,plot_f_1,rgb,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,mser_min_range,mser_max_range,close_enough );
        lmm = [lmm ;ri.PixelList ];
        

end
fprintf('computing for H \n');
[MserRegions_H,~] = detectMSERFeatures(H, 'RegionAreaRange',[mser_min_range ,mser_max_range],'ThresholdDelta',ThresholdDelta);
[~, rH] =  filter_doublon_intra(MserRegions_H,plot_f_1,rgb,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,mser_min_range,mser_max_range,close_enough );
lmm = [lmm ;rH.PixelList ];

fprintf('computing for E \n');
[MserRegions_E,~] = detectMSERFeatures(E, 'RegionAreaRange',[mser_min_range ,mser_max_range],'ThresholdDelta',ThresholdDelta);
[~, rE] =  filter_doublon_intra(MserRegions_E,plot_f_1,rgb,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,mser_min_range,mser_max_range,close_enough );
lmm = [lmm ;rE.PixelList ];


concf_MSER_after_inter = MSERRegions(lmm);
regions_after_inter_filter = filter_doublon_inter(concf_MSER_after_inter ,plot_f_1,rgb,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,mser_min_range,mser_max_range,close_enough );
regions_after_inter_filter = regions_after_inter_filter(:,1);

regions_candidates_after_hne_filter_pixel_list = candidates_mser_regions_after_hne_filter(H, E,MSERRegions(regions_after_inter_filter), HnE_quantile_treshold, []);
[s1,~]= size(regions_candidates_after_hne_filter_pixel_list); 
if s1 ~= 0
    
    
    
    regions_candidates_after_hne_filter = MSERRegions(regions_candidates_after_hne_filter_pixel_list);      
    figure('Name','cell segmentation vizualisation');imshow(rgb); hold on;
    hp = impixelinfo;
    plot(regions_candidates_after_hne_filter,'showPixelList',false,'showEllipses',true);
    hold off;


    %post process 

    fprintf('postprocess \n');
    [s1,~] = size(regions_candidates_after_hne_filter_pixel_list);
    [bg_row,bg_col,~] = size(rgb);
        New_regions_candidates_after_hne_filter_pixel_list = [];

    for j=1:s1
        cj  =regions_candidates_after_hne_filter_pixel_list{j};
        v = zeros(bg_row, bg_col);
        [l1, ~] = size(cj);
        for i=1:l1
            v(cj(i,1),cj(i,2)) = true;
        end
        t = bwmorph(v, 'thicken');
        c = imfill(t ,'holes');
        [c1,c2] = size(c);
        ci = [];
        for k =1:c1
            for l =1:c2
                if c(k,l)
                    ci = [ ci ; k,l];
                end
            end
        end
        New_regions_candidates_after_hne_filter_pixel_list = [New_regions_candidates_after_hne_filter_pixel_list; {ci}];
    end
    
    pred_list = New_regions_candidates_after_hne_filter_pixel_list;
else 
    pred_list = [];
end

end
%}
