%
%x20 :  100 , 4000
%x40 : 200 , 8000


% tissue_area_aray_of_candidate_cells = [];

%###############    START      ######################
    %1 hyper para 
%Set initial default parameter
%default cut in spectrum dimension:
bottom_wl_treshold = 4;
top_wl_treshold = 30;
%default filter smoothing: medfilt2. todo: choose which default
%default nbr of clusters in weighted kmean:
nbr_of_cluster = 20;

%TODOOOOO  generlaize relevant cluster by tissue



%##############for tissue: 'x40_r1_r2' of C161VA_x20_X700Y600_
%{
imaging = 'x40';
fprintf('x40_r1_r2 \n');
tissue_name ='x40_r1_r2';
weighted_kmean_name = 'x24wkm_20_kmeans_and_D.mat';

crop_area_1 = [500 1000  (1000-500) (1500-1000)];
crop_area_2 = [1750 580  (2150-1750) (820-580)];
crop_area_3 = [ 1500 1000  (2000-1500) (1500-1000)];
crop_area_4 = [ 3400 1000  (3800-3400) (1400-1000)];
all = [];

seleted_area_list = [ crop_area_1 ; crop_area_2 ; crop_area_3 ; crop_area_4 ];
selected_area_name_list  = [ 'area_1';  'area_2';  'area_3'; 'area_4'];

%selected_area = seleted_area_list(4,:);
%selected_area = all;


relevant_clusters = [0.21, 0.05, 0.37];  



%###########  select area we work on - if all the tissue m select area = all
%folder_path  = '/Users/lior/Desktop/Image & Analysis /cells_label/'+ string(tissue_name)+string('/')+selected_area_name+'/';


[data , lambda_mat , wl_mat ,rgb,  distance_map_list] = initialization_for_one_tissue(tissue_name, weighted_kmean_name, bottom_wl_treshold, top_wl_treshold, nbr_of_cluster, relevant_clusters); 



%HNE + agadire : new_tissue_area_aray_of_candidate_cells



%############## get HnE
[full_H, full_E] = get_HnE(wl_mat,lambda_mat);


%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------


%%%%% STEP AGADIR -  pr chaque area de tissue,  recupere les candidat_list

if imaging == 'x40'
    mser_min_range = 200;
    mser_max_range = 8000;
else
    mser_min_range = 100;
    mser_max_range = 4000; 
end
nbr_of_areas_by_tissue = 4;
%tissue_name_list = ['c161'; 'x40_r1_r2'];
tissue_name_list = [tissue_name];
[nbr_of_tissues, ~] = size(tissue_name_list);
%tissue_area_aray_of_candidate_cells = [];
%tissue_area_aray_of_candidate_cells = load('tissue_area_aray_of_candidate_cells.mat').tissue_area_aray_of_candidate_cells;

for t= 1:nbr_of_tissues
    %here insert all tissue calc including hnE
    for a=1: nbr_of_areas_by_tissue
        fprintf('area %d \n ' , a);
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
        min_range = mser_min_range;
        max_range = mser_max_range;
        if imaging == 'x40'
            close_enough = 40;
        else 
            close_enough = 20; 
        end

        [MserRegions_1,~] = detectMSERFeatures(distance_map_list_a1, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r1] =  filter_doublon_intra(MserRegions_1,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_2,~] = detectMSERFeatures(distance_map_list_a2, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r2] =  filter_doublon_intra(MserRegions_2,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_3,~] = detectMSERFeatures(distance_map_list_a3, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r3] =  filter_doublon_intra(MserRegions_3,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_H,~] = detectMSERFeatures(H, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, rH] =  filter_doublon_intra(MserRegions_H,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_E,~] = detectMSERFeatures(E, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, rE] =  filter_doublon_intra(MserRegions_E,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );


        concf_MSER_after_inter = MSERRegions([r1.PixelList;r2.PixelList;r3.PixelList; rH.PixelList; rE.PixelList]);


        %##############  Filter_1 : Doublons: removing closes/included ellipses

        [~,regions_after_inter_filter] = filter_doublon_inter(concf_MSER_after_inter ,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );
        %~~~~~~~~ ie F
        
        %vizualise:
        %{
        figure('Name',tissue_name + " - " + selected_area_name);imshow(background); hold on;
        hp = impixelinfo;
        plot(regions_after_inter_filter,'showPixelList',false,'showEllipses',true);
        hold off;
        %}
        
        tissue_area_aray_of_candidate_cells = [tissue_area_aray_of_candidate_cells ; {regions_after_inter_filter}, {full_H}, {full_E}, {folder_path}, {selected_area}  ] ;
    end
end
%save('old_tissue_area_aray_of_candidate_cells', 'tissue_area_aray_of_candidate_cells');
%save('tissue_area_aray_of_candidate_cells', 'tissue_area_aray_of_candidate_cells');
%tissue_area_aray_of_candidate_cells = load('tissue_area_aray_of_candidate_cells.mat').tissue_area_aray_of_candidate_cells;

%}
















%#####################for tissue '3.rawstack'
%{
fprintf( "3.rawstack \n ");
tissue_name = "3.rawstack";
struct_path = fullfile(pwd, tissue_name + ".mat");
imaging = 'x20';
crop_area_1 = [600 300 200 300];
crop_area_2 = [400 800 287 146];
crop_area_3 = [1000 900 222 100];
crop_area_4 = [700 500 100 100];
all = [];

seleted_area_list = [ crop_area_1 ; crop_area_2 ; crop_area_3 ; crop_area_4 ];
selected_area_name_list  = [ 'area_1';  'area_2';  'area_3'; 'area_4'];

%if smooth cut wl and wkmean not computed already:
%{
data = load(struct_path);
lambda_mat = data.lambda;
wl_mat = data.spec;
RGB = data.RGB;
%imtool(RGB);
%}



        

        
%-----------------------------------------------------------------

%##############  Cutting and smothing  wl from spec
 %cut smooth wl + kmeans : compute or load       
 
%TOdO : filter out of wl_cut + generic filter instead of medfilt2

%  cut and same smoothing to all WL , and  save ?  the new spectrum
%{
smoothed_and_wl_shaped = medfilt2(read_unique_wl(wl_mat , lambda_mat, bottom_wl_treshold));
for i = bottom_wl_treshold+1 :top_wl_treshold
single_wl = read_unique_wl(wl_mat , lambda_mat, i);
smoothed_wl= medfilt2(single_wl);
smoothed_and_wl_shaped = cat(3,smoothed_and_wl_shaped, smoothed_wl);

end
wl_mat = smoothed_and_wl_shaped;
lambda_mat = lambda_mat(bottom_wl_treshold:top_wl_treshold);
spec = wl_mat;
lambda = lambda_mat;
size(smoothed_and_wl_shaped);
%VISU
cutted_and_smoothed_spectral_img =  spectoRGB(smoothed_and_wl_shaped,lambda_mat, 0,'name');
cutted_and_smoothed_spectral_img_formatted =  uint8(255 * mat2gray(cutted_and_smoothed_spectral_img));
imshow(cutted_and_smoothed_spectral_img_formatted)



%SAVING: 
%inutile%smoothed_and_wl_shaped  = cast(new_mat, 'double');
name = string(tissue_name )+ "_smoothed_and_wl_shaped.mat" ;
save(name,'spec', 'lambda', 'RGB');

%}
%or load cutted and smoothed wl 
data = load(fullfile(pwd,tissue_name+"_smoothed_and_wl_shaped.mat" ));
rgb = data.RGB;
wl_mat = data.spec;
lambda_mat = data.lambda;



%##############  II kmeans

%TODO: full path as new last param in wkm fct

%compute weighted kmeans
%{
input_spectral_img = wl_mat; %wl_mat
%input_spectral_img = imcrop(input_spectral_img,[2761 686  (3781-2761) (1409-686)]);
number_of_clusters = nbr_of_cluster;

save_kmean_img=1; %1 if want to save img
name_of_input_img = tissue_name+"_smoothed_and_wl_shaped";

[W_kmean_mat_output, D]=weighted_kmeans_algo(input_spectral_img,number_of_clusters, save_kmean_img, name_of_input_img,pwd);

%}

%or load weighted kmeans 

wkm= load('3.rawstack_smoothed_and_wl_shaped_20_kmeans_and_D.mat');
W_kmean_mat_output = wkm.kmean_mat_output;
D = wkm.D;


%imtool(wkm.kmean_mat_output);
%km = uint8(255*(wkm.kmean_mat_output));

%-----------------------------------------------------------------

%##############  III select good fit cluster 

%ask user to select relevant cluster values
% when watching %imtool(wkm.kmean_mat_output);
%imtool(W_kmean_mat_output);


relevant_cluster = [.05 , .26, .84 ]; % for tissue 3.rawstack ie '3.rawstack_smoothed_and_wl_shaped_20_kmeans_and_D.mat'

%relevant_cluster conversion : from value([0:1] bcs in graylvl) to dim (nbr_of_cluster)
[~ , s2]= size(relevant_cluster);
dim_elt_in_km = [];
for elt = 1:s2
    dim_elt_in_km = [dim_elt_in_km, round(nbr_of_cluster * relevant_cluster(elt))];
end


%visualize selected relevant clusters
%{
[a,b]=size(W_kmean_mat_output);
F = reshape(D, [a,b,number_of_clusters]);
for i=1:s2
    imtool(mat2gray(F(:,:,dim_elt_in_km(i)))); %dbg here 
end
%}

%----------------------------------------------------------------------
%##############  IV DISTANCE MATRIX FROM KMEANS

%Load, preprocess, inverse & normalize, distance of each pixel from
%centroids and get the relevant clusters.

%inverse and normalise:
X =(1./D);
XN = normalize(X,2,'range');

%reshape/preprocess
[a,b]=size(W_kmean_mat_output);
distance_map_list= zeros(a,b, s2);
for elt= 1:s2  
distance_map_list(:,:,elt) = reshape(XN(:,dim_elt_in_km(elt)), [a,b]);
end

%----------------------------------------------------------------------
%mser_min_range = 200;
%mser_max_range = 5000;
%----------------------------------------------------------------------

%HNE + agadire : new_tissue_area_aray_of_candidate_cells

%############## get HnE
% if not working , le mettre plus bas.
[full_H, full_E] = get_HnE(wl_mat,lambda_mat);


%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------


%%%%% STEP AGADIR -  pr chaque area de tissue,  recupere les candidat_list
if imaging == 'x40'
    mser_min_range = 200;
    mser_max_range = 8000;
else
    mser_min_range = 100;
    mser_max_range = 4000; 
end
nbr_of_areas_by_tissue = 4;
tissue_name_list = [tissue_name];
[nbr_of_tissues, ~] = size(tissue_name_list);
%tissue_area_aray_of_candidate_cells = [];
%tissue_area_aray_of_candidate_cells = load('tissue_area_aray_of_candidate_cells.mat').tissue_area_aray_of_candidate_cells;

for t= 1:nbr_of_tissues
    for a=1: nbr_of_areas_by_tissue
        fprintf('area %d \n ' , a);
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
        min_range = mser_min_range;
        max_range = mser_max_range;
        if imaging == 'x40'
            close_enough = 40;
        else 
            close_enough = 20; 
        end

        [MserRegions_1,~] = detectMSERFeatures(distance_map_list_a1, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r1] =  filter_doublon_intra(MserRegions_1,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_2,~] = detectMSERFeatures(distance_map_list_a2, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r2] =  filter_doublon_intra(MserRegions_2,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_3,~] = detectMSERFeatures(distance_map_list_a3, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r3] =  filter_doublon_intra(MserRegions_3,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );


        [MserRegions_H,~] = detectMSERFeatures(H, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, rH] =  filter_doublon_intra(MserRegions_H,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_E,~] = detectMSERFeatures(E, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, rE] =  filter_doublon_intra(MserRegions_E,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );


        concf_MSER_after_inter = MSERRegions([r1.PixelList;r2.PixelList;r3.PixelList; rH.PixelList; rE.PixelList]);


        %##############  Filter_1 : Doublons: removing closes/included ellipses

        [~,regions_after_inter_filter] = filter_doublon_inter(concf_MSER_after_inter ,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );
        %~~~~~~~~ ie F
        
        %vizualise:
        %{
        figure('Name',tissue_name + " - " + selected_area_name);imshow(background); hold on;
        hp = impixelinfo;
        plot(regions_after_inter_filter,'showPixelList',false,'showEllipses',true);
        hold off;
        %}
        
        tissue_area_aray_of_candidate_cells = [tissue_area_aray_of_candidate_cells ; {regions_after_inter_filter}, {full_H}, {full_E}, {folder_path}, {selected_area}  ] ;
    end
end
%save('old_tissue_area_aray_of_candidate_cells', 'tissue_area_aray_of_candidate_cells');
%save('tissue_area_aray_of_candidate_cells', 'tissue_area_aray_of_candidate_cells');
%tissue_area_aray_of_candidate_cells = load('tissue_area_aray_of_candidate_cells.mat').tissue_area_aray_of_candidate_cells;

%}











%#####################for tissue 'C083HA'
%{
%previous tissue "C083HA_0_1800_0"
fprintf('C083HA \n');
imaging = 'x40';

tissue_name = 'C083HA';
struct_path = 'C083HA_0_1800_0_clean.mat';



crop_area_1 = [560 1 440 500];
crop_area_2 = [4300 150 300 350];
crop_area_3 = [5800 300 300 150];
crop_area_4 = [3850 850 600 150];
all = [];

seleted_area_list = [ crop_area_1 ; crop_area_2 ; crop_area_3 ; crop_area_4 ];
selected_area_name_list  = [ 'area_1';  'area_2';  'area_3'; 'area_4'];




%##############  Cutting and smothing  wl from spec
 %cut smooth wl + kmeans : compute or load       
 
%TOdO : filter out of wl_cut + generic filter instead of medfilt2

%  cut and same smoothing to all WL , and  save ?  the new spectrum
%{
smoothed_and_wl_shaped = medfilt2(read_unique_wl(wl_mat , lambda_mat, bottom_wl_treshold));
for i = bottom_wl_treshold+1 :top_wl_treshold
single_wl = read_unique_wl(wl_mat , lambda_mat, i);
smoothed_wl= medfilt2(single_wl);
smoothed_and_wl_shaped = cat(3,smoothed_and_wl_shaped, smoothed_wl);

end
wl_mat = smoothed_and_wl_shaped;
lambda_mat = lambda_mat(bottom_wl_treshold:top_wl_treshold);
spec = wl_mat;
lambda = lambda_mat;
size(smoothed_and_wl_shaped);
%VISU
cutted_and_smoothed_spectral_img =  spectoRGB(smoothed_and_wl_shaped,lambda_mat, 0,'name');
cutted_and_smoothed_spectral_img_formatted =  uint8(255 * mat2gray(cutted_and_smoothed_spectral_img));
imshow(cutted_and_smoothed_spectral_img_formatted)



%SAVING: 
%inutile%smoothed_and_wl_shaped  = cast(new_mat, 'double');
%name = string(tissue_name )+ "_smoothed_and_wl_shaped.mat" ;
%save(name,'spec', 'lambda', 'RGB');

%}
%or load cutted and smoothed wl 
data = load(struct_path);
lambda_mat = data.lambda;
wl_mat = data.spec;
rgb = data.RGB;
%imtool(rgb);


%##############  II kmeans

%TODO: full path as new last param in wkm fct
%{
%compute weighted kmeans

input_spectral_img = wl_mat; %wl_mat
%input_spectral_img = imcrop(input_spectral_img,[2761 686  (3781-2761) (1409-686)]);
number_of_clusters = nbr_of_cluster;

save_kmean_img=1; %1 if want to save img
name_of_input_img = tissue_name+"_smoothed_and_wl_shaped";

[W_kmean_mat_output, D]=weighted_kmeans_algo(input_spectral_img,number_of_clusters, save_kmean_img, name_of_input_img,pwd);

%}

%or load weighted kmeans 
wkm = load('C083HA_0_1800_0_clean_20_kmeans_and_D.mat');
W_kmean_mat_output = wkm.kmean_mat_output;
D = wkm.D;



%imtool(wkm.kmean_mat_output);
%km = uint8(255*(wkm.kmean_mat_output));

%-----------------------------------------------------------------

%##############  III select good fit cluster 

%ask user to select relevant cluster values
% when watching %imtool(wkm.kmean_mat_output);
%imtool(W_kmean_mat_output);

relevant_cluster = [.68 , .95 ]; % for tissue C083HA_0_1800_0_clean.mat

%relevant_cluster conversion : from value([0:1] bcs in graylvl) to dim (nbr_of_cluster)
[~ , s2]= size(relevant_cluster);
dim_elt_in_km = [];
for elt = 1:s2
    dim_elt_in_km = [dim_elt_in_km, round(nbr_of_cluster * relevant_cluster(elt))];
end


%visualize selected relevant clusters
%{
[a,b]=size(W_kmean_mat_output);
F = reshape(D, [a,b,number_of_clusters]);
for i=1:s2
    imtool(mat2gray(F(:,:,dim_elt_in_km(i)))); %dbg here 
end
%}

%----------------------------------------------------------------------
%##############  IV DISTANCE MATRIX FROM KMEANS

%Load, preprocess, inverse & normalize, distance of each pixel from
%centroids and get the relevant clusters.

%inverse and normalise:
X =(1./D);
XN = normalize(X,2,'range');

%reshape/preprocess
[a,b]=size(W_kmean_mat_output);
distance_map_list= zeros(a,b, s2);
for elt= 1:s2  
distance_map_list(:,:,elt) = reshape(XN(:,dim_elt_in_km(elt)), [a,b]);
end

%----------------------------------------------------------------------
%mser_min_range = 250;
%mser_max_range = 10000;
%----------------------------------------------------------------------


%HNE + agadire : new_tissue_area_aray_of_candidate_cells

%############## get HnE
% if not working , le mettre plus bas.
[full_H, full_E] = get_HnE(wl_mat,lambda_mat);


%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------


%%%%% STEP AGADIR -  pr chaque area de tissue,  recupere les candidat_list
if imaging == 'x40'
    mser_min_range = 200;
    mser_max_range = 8000;
else
    mser_min_range = 100;
    mser_max_range = 4000; 
end

nbr_of_areas_by_tissue = 4;
%tissue_name_list = ['c161'; 'x40_r1_r2'];
tissue_name_list = [tissue_name];
[nbr_of_tissues, ~] = size(tissue_name_list);
%tissue_area_aray_of_candidate_cells = [];
%tissue_area_aray_of_candidate_cells = load('tissue_area_aray_of_candidate_cells.mat').tissue_area_aray_of_candidate_cells;

for t= 1:nbr_of_tissues
    %here insert all tissue calc including hnE
    for a=1: nbr_of_areas_by_tissue
        fprintf('area %d \n ' , a);
        selected_area = seleted_area_list(a,:);
        selected_area_name =  selected_area_name_list(a,:);
        folder_path = pwd  + "/cells_label/" + string(tissue_name) + "/" +selected_area_name+ "/" ;

        [~ , s2_sa] = size(selected_area);
        if s2_sa == 0
            background = rgb;
            distance_map_list_a1 = distance_map_list(:,:,1);
            distance_map_list_a2 = distance_map_list(:,:,2);
            H = full_H;
            E = full_E;
        else
            %distance_map_list_full_a1 = distance_map_list(:,:,1);
            background = imcrop(rgb, selected_area);
            distance_map_list_a1 = imcrop(distance_map_list(:,:,1), selected_area);
            distance_map_list_a2 = imcrop(distance_map_list(:,:,2), selected_area);

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
        min_range = mser_min_range;
        max_range = mser_max_range;
        if imaging == 'x40'
            close_enough = 40;
        else 
            close_enough = 20; 
        end

        [MserRegions_1,~] = detectMSERFeatures(distance_map_list_a1, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r1] =  filter_doublon_intra(MserRegions_1,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_2,~] = detectMSERFeatures(distance_map_list_a2, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r2] =  filter_doublon_intra(MserRegions_2,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );


        [MserRegions_H,~] = detectMSERFeatures(H, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, rH] =  filter_doublon_intra(MserRegions_H,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_E,~] = detectMSERFeatures(E, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, rE] =  filter_doublon_intra(MserRegions_E,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );


        concf_MSER_after_inter = MSERRegions([r1.PixelList;r2.PixelList; rH.PixelList; rE.PixelList]);


        %##############  Filter_1 : Doublons: removing closes/included ellipses

        [~,regions_after_inter_filter] = filter_doublon_inter(concf_MSER_after_inter ,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );
        %~~~~~~~~ ie F
        
        %vizualise:
        %{
        figure('Name',tissue_name + " - " + selected_area_name);imshow(background); hold on;
        hp = impixelinfo;
        plot(regions_after_inter_filter,'showPixelList',false,'showEllipses',true);
        hold off;
        %}
        
        tissue_area_aray_of_candidate_cells = [tissue_area_aray_of_candidate_cells ; {regions_after_inter_filter}, {full_H}, {full_E}, {folder_path}, {selected_area}  ] ;
    end
end
%save('old_tissue_area_aray_of_candidate_cells', 'tissue_area_aray_of_candidate_cells');
%save('tissue_area_aray_of_candidate_cells', 'tissue_area_aray_of_candidate_cells');
%tissue_area_aray_of_candidate_cells = load('tissue_area_aray_of_candidate_cells.mat').tissue_area_aray_of_candidate_cells;

%}








%dbln x40
%{
%###################### for tissue C161VA_x20_X700Y600
%{
fprintf('C161VA_x20_X700Y600 \n');
imaging  = 'x20';
tissue_name = 'C161VA_x20_X700Y600';
struct_path = 'C161VA_x20_X700Y600.mat';



crop_area_1 = [1700 700 300 300];
crop_area_2 = [1200 1500 300 200];
crop_area_3 = [500 800 150 300];
crop_area_4 = [800 550 400 150];
all = [];

seleted_area_list = [ crop_area_1 ; crop_area_2 ; crop_area_3 ; crop_area_4 ];
selected_area_name_list  = [ 'area_1';  'area_2';  'area_3'; 'area_4'];





%##############  Cutting and smothing  wl from spec
 %cut smooth wl + kmeans : compute or load       
 
%TOdO : filter out of wl_cut + generic filter instead of medfilt2

%  cut and same smoothing to all WL , and  save ?  the new spectrum
%{
smoothed_and_wl_shaped = medfilt2(read_unique_wl(wl_mat , lambda_mat, bottom_wl_treshold));
for i = bottom_wl_treshold+1 :top_wl_treshold
single_wl = read_unique_wl(wl_mat , lambda_mat, i);
smoothed_wl= medfilt2(single_wl);
smoothed_and_wl_shaped = cat(3,smoothed_and_wl_shaped, smoothed_wl);
end
wl_mat = smoothed_and_wl_shaped;
lambda_mat = lambda_mat(bottom_wl_treshold:top_wl_treshold);
spec = wl_mat;
lambda = lambda_mat;
size(smoothed_and_wl_shaped);

%Vizualisation
%cutted_and_smoothed_spectral_img =  spectoRGB(smoothed_and_wl_shaped,lambda_mat, 0,'name');
%cutted_and_smoothed_spectral_img_formatted =  uint8(255 * mat2gray(cutted_and_smoothed_spectral_img));
%imshow(cutted_and_smoothed_spectral_img_formatted);

%SAVING: 
%inutile%smoothed_and_wl_shaped  = cast(new_mat, 'double');
%name = string(tissue_name )+ "_smoothed_and_wl_shaped.mat" ;
%save(name,'spec', 'lambda', 'RGB');
%}

%or load cutted and smoothed wl 

data = load('C161VA_x20_X700Y600_smoothed_and_wl_shaped.mat');
lambda_mat = data.lambda;
wl_mat = data.spec;
rgb = data.RGB;





%##############  II kmeans

%TODO: full path as new last param in wkm fct

%compute weighted kmeans
%{
input_spectral_img = wl_mat; %wl_mat
%input_spectral_img = imcrop(input_spectral_img,[2761 686  (3781-2761) (1409-686)]);
number_of_clusters = nbr_of_cluster;

save_kmean_img=1; %1 if want to save img
name_of_input_img = tissue_name+"_smoothed_and_wl_shaped";

[W_kmean_mat_output, D]=weighted_kmeans_algo(input_spectral_img,number_of_clusters, save_kmean_img, name_of_input_img,pwd);
%}

%or load weighted kmeans 

wkm = load('C161VA_x20_X700Y600_smoothed_and_wl_shaped_20_kmeans_and_D.mat');
W_kmean_mat_output = wkm.kmean_mat_output;
D = wkm.D;



%imtool(wkm.kmean_mat_output);
%km = uint8(255*(wkm.kmean_mat_output));

%-----------------------------------------------------------------

%##############  III select good fit cluster 

%ask user to select relevant cluster values
% when watching %imtool(wkm.kmean_mat_output);
%imtool(W_kmean_mat_output);


relevant_cluster = [.89 ,.16 ,.05 , 0.32 ]; % for tissue 'C161VA_x20_X700Y600'.mat

%relevant_cluster conversion : from value([0:1] bcs in graylvl) to dim (nbr_of_cluster)
[~ , s2]= size(relevant_cluster);
dim_elt_in_km = [];
for elt = 1:s2
    dim_elt_in_km = [dim_elt_in_km, round(nbr_of_cluster * relevant_cluster(elt))];
end


%visualize selected relevant clusters
%{
[a,b]=size(W_kmean_mat_output);
F = reshape(D, [a,b,number_of_clusters]);
for i=1:s2
    imtool(mat2gray(F(:,:,dim_elt_in_km(i)))); %dbg here 
end
%}

%----------------------------------------------------------------------
%##############  IV DISTANCE MATRIX FROM KMEANS

%Load, preprocess, inverse & normalize, distance of each pixel from
%centroids and get the relevant clusters.

%inverse and normalise:
X =(1./D);
XN = normalize(X,2,'range');

%reshape/preprocess
[a,b]=size(W_kmean_mat_output);
distance_map_list= zeros(a,b, s2);
for elt= 1:s2  
distance_map_list(:,:,elt) = reshape(XN(:,dim_elt_in_km(elt)), [a,b]);
end

%----------------------------------------------------------------------
mser_min_range = 100;
mser_max_range = 1500;
%----------------------------------------------------------------------
%}
%HNE + agadir : tissue_area_aray_of_candidate_cells
%{
%############## get HnE
[full_H, full_E] = get_HnE(wl_mat,lambda_mat);


%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------


%%%%% STEP AGADIR -  pr chaque area de tissue,  recupere les candidat_list
nbr_of_areas_by_tissue = 4;
%tissue_name_list = ['c161'; 'x40_r1_r2'];
tissue_name_list = [tissue_name];
[nbr_of_tissues, ~] = size(tissue_name_list);
%tissue_area_aray_of_candidate_cells = [];
%tissue_area_aray_of_candidate_cells = load('tissue_area_aray_of_candidate_cells.mat').tissue_area_aray_of_candidate_cells;

for t= 1:nbr_of_tissues
    for a=1: nbr_of_areas_by_tissue
        fprintf('area %d \n ' , a);
        selected_area = seleted_area_list(a,:);
        selected_area_name =  selected_area_name_list(a,:);
        folder_path = pwd  + "/cells_label/" + string(tissue_name) + "/" +selected_area_name+ "/" ;

        [~ , s2_sa] = size(selected_area);
        if s2_sa == 0
            background = rgb;
            distance_map_list_a1 = distance_map_list(:,:,1);
            distance_map_list_a2 = distance_map_list(:,:,2);
            distance_map_list_a3 = distance_map_list(:,:,3);
            distance_map_list_a4 = distance_map_list(:,:,4);
            H = full_H;
            E = full_E;
        else
            %distance_map_list_full_a1 = distance_map_list(:,:,1);
            background = imcrop(rgb, selected_area);
            distance_map_list_a1 = imcrop(distance_map_list(:,:,1), selected_area);
            distance_map_list_a2 = imcrop(distance_map_list(:,:,2), selected_area);
            distance_map_list_a3 = imcrop(distance_map_list(:,:,3), selected_area);
            distance_map_list_a4 = imcrop(distance_map_list(:,:,4), selected_area);
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
        min_range = mser_min_range;
        max_range = mser_max_range;
        if imaging == 'x40'
            close_enough = 40;
        else 
            close_enough = 20; 
        end

        [MserRegions_1,~] = detectMSERFeatures(distance_map_list_a1, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r1] =  filter_doublon_intra(MserRegions_1,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_2,~] = detectMSERFeatures(distance_map_list_a2, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r2] =  filter_doublon_intra(MserRegions_2,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_3,~] = detectMSERFeatures(distance_map_list_a3, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r3] =  filter_doublon_intra(MserRegions_3,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_4,~] = detectMSERFeatures(distance_map_list_a4, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, r4] =  filter_doublon_intra(MserRegions_4,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_H,~] = detectMSERFeatures(H, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, rH] =  filter_doublon_intra(MserRegions_H,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );

        [MserRegions_E,~] = detectMSERFeatures(E, 'RegionAreaRange',[min_range ,max_range],'ThresholdDelta',ThresholdDelta);
        [~, rE] =  filter_doublon_intra(MserRegions_E,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );


        concf_MSER_after_inter = MSERRegions([r1.PixelList;r2.PixelList;r3.PixelList; r4.PixelList; rH.PixelList; rE.PixelList]);


        %##############  Filter_1 : Doublons: removing closes/included ellipses

        [~,regions_after_inter_filter] = filter_doublon_inter(concf_MSER_after_inter ,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );
        %~~~~~~~~ ie F
        
        %vizualise:
        
        figure('Name',tissue_name + " - " + selected_area_name);imshow(background); hold on;
        hp = impixelinfo;
        plot(regions_after_inter_filter,'showPixelList',false,'showEllipses',true);
        hold off;
        
        
        tissue_area_aray_of_candidate_cells = [tissue_area_aray_of_candidate_cells ; {regions_after_inter_filter}, {full_H}, {full_E}, {folder_path}, {selected_area}  ] ;
    end
end
%save('new_tissue_area_aray_of_candidate_cells' ,'new_tissue_area_aray_of_candidate_cells');
%save('old_tissue_area_aray_of_candidate_cells', 'tissue_area_aray_of_candidate_cells');
%save('tissue_area_aray_of_candidate_cells', 'tissue_area_aray_of_candidate_cells');
%tissue_area_aray_of_candidate_cells = load('tissue_area_aray_of_candidate_cells.mat').tissue_area_aray_of_candidate_cells;

%}
%}

%}
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%tissue_area_aray_of_candidate_cells = load('new_tissue_area_aray_of_candidate_cells.mat').new_tissue_area_aray_of_candidate_cells;

%%%%%%% STEP BARUCH : get confumat 

[q,~] = size(tissue_area_aray_of_candidate_cells);
%q = 1*4 + 1*4;
confumat = zeros(q, 21, 21,3);
for iou_treshold= 0:.05:1
    fprintf('iou_treshold %d \n', iou_treshold);
    for HnE_quantile_treshold=0:.05:1
        for zone = 1:q
            candidates_regions = tissue_area_aray_of_candidate_cells{zone, 1};
          %{  
            for region_index=1:candidates_regions.Count
                regions_pixels = candidates_regions.PixelList(region_index);
                if max(max(regions_pixels)) > 301 && zone==13
                    lkd=0;
                end
            end
            %}
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
title('Precision-Recall graph for different IoU threshold');
hold on;
[ss,~] = size(precision_by_iou_th_list);
for i=2:ss
    plot(recall_by_iou_th_list(i,:) , precision_by_iou_th_list(i,:))
end
hold off; legend(lg, 'Location','southwest');

%-----------------------------------------------------------------
%}
