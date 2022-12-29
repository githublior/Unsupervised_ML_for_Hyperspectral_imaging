%PIPE
%{
- option whitening  + select  white rectangles
- preprocess :  remove some wl + filter . which  ?
- wkm
- select cluster - (automatise conversion)
- filtered selection to vizualize and validate our cluster selection
- mser
-filter 1

%}
%TOdo 
%{
- im filter sortir du wk coupage et generaliser le type de filter
 - whitening apres coupe
- on veut du high recall avant tout . jouer avec fpr tpr pour decider best
params.
- requierd un exemple taggé




%}
%start /vizu
%{
 %prompt= 'input struct .mat full path (contains at least .spec and .lambda matrixes)';
 %struct_path=input(prompt);
 %wl_struct = load(struct_path);


 %if create rgb: save_SpecToRgb = true:
 %prompt = 'spec to rgb tissue name?';
 %name = input(prompt)
 %save_SpecToRgb = spectoRgb(wl_struct.spec,wl_struct.lambda, save_SpecToRgb,name);
 %imshow(save_SpecToRgb)
% hf. ... ?
%if whitening. ....

%}



%-----------------------------------------------------------------
%TODO : whitening ? preprocess: whitening cyclic noise ?
%see start of b5
%I preprocess

%Set initial default parameter

%default cut in spectrum dimension:
bottom_wl_treshold = 4;
top_wl_treshold = 30;
%default filter smoothing: medfilt2. todo: choose which default


%default nbr of clusters in weighted kmean:
nbr_of_cluster = 20;

%-----------------------------------------------------------------


%############## I initialization -uncoment

%previous tissue x40_r1_r2
%{
% load tissue we work on 

%wetransfer_3-rawstack-jpg_2021-08-30_1417/
tissue_name ='x40_r1_r2';


%TO CHECK NEXT LINE 
struct_path = pwd  + string('/') + string(tissue_name) + '.mat';
%struct_path = '/Users/lior/Desktop/Image & Analysis /' +  string(tissue_name) + '.mat';
data = load(struct_path);

%imtool(data.rgb)
lambda_mat = data.lambda;
rgb = data.rgb;
wl_mat = data.spec;
[a,b, ~ ] = size(wl_mat);


%-----------------------------------------------------------------

%############## get HnE


% if not working , le mettre plus bas.

[full_H, full_E] = get_HnE(wl_mat,lambda_mat);
%}

%previous tissue "3.rawstack"
%{
tissue_name = "3.rawstack";
struct_path = fullfile(pwd, tissue_name + ".mat");

crop_area_1 = [600 300 200 300];
crop_area_2 = [400 800 287 146];
crop_area_3 = [1000 900 222 100];
crop_area_4 = [700 500 100 100];

selected_area = crop_area_1;
selected_area_name = 'area_1';

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


% relevant_cluster = [0.21, 0.05, 0.37]; % for tissue x40_r1_r2
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
%}
%----------------------------------------------------------------------
%----------------------------------------------------------------------




%previous tissue "C083HA_0_1800_0"
%{

struct_path = 'C083HA_0_1800_0_clean.mat';



crop_area_1 = [560 1 440 500];
crop_area_2 = [4300 150 300 350];
crop_area_3 = [5800 300 300 150];
crop_area_4 = [3850 850 600 150];





%imtool(RGB);







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
RGB = data.RGB;





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


% relevant_cluster = [0.21, 0.05, 0.37]; % for tissue x40_r1_r2
% relevant_cluster = [.05 , .26, .84 ]; % for tissue 3.rawstack ie '3.rawstack_smoothed_and_wl_shaped_20_kmeans_and_D.mat'
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
%----------------------------------------------------------------------


%}




%previous tissue "C161VA_x20_X700Y600"
%{
tissue_name = 'C161VA_x20_X700Y600';
struct_path = 'C161VA_x20_X700Y600.mat';
data = load(struct_path);

%imtool(data.rgb)
lambda_mat = data.lambda;
RGB = data.RGB;
wl_mat = data.spec;
%imtool(RGB);


crop_area_1 = [1700 700 300 300];
crop_area_2 = [1200 1500 300 200];
crop_area_3 = [500 800 150 300];
crop_area_4 = [800 550 400 150];




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
RGB = data.RGB;





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
%----------------------------------------------------------------------


%}



%HnE
%[full_H, full_E] = get_HnE(wl_mat,lambda_mat);
%----------------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%-----------------------------------------------------------------




%###############    after fine tuning , use parameter to show results      ######################


%-----------------------------------------------------------------

%##############  select area we work on - if all the tissue m select area = all
%{

crop_area_1 = [500 1000  (1000-500) (1500-1000)];
crop_area_2 = [1750 580  (2150-1750) (820-580)];
crop_area_3 = [ 1500 1000  (2000-1500) (1500-1000)];
crop_area_4 = [ 3400 1000  (3800-3400) (1400-1000)];
all = [];
%crop_area_11 = [500+370 1000+400  (1000-(500+370)) (1500-(1000+400))];

selected_area =crop_area_4;
selected_area_name = 'area_4';

%selected_area = all;

     
%folder_path  = '/Users/lior/Desktop/Image & Analysis /cells_label/'+ string(tissue_name)+string('/')+selected_area_name+'/';

folder_path = pwd  + "/cells_label/" + string(tissue_name) + "/" +selected_area_name+ "/" ;

%}
%-----------------------------------------------------------------


  %      !! regions_candidates_after_hne_filter_pixel_list = candidates_mser_regions_after_hne_filter(full_H, full_E,regions_after_inter_filter, HnE_quantile_treshold, selected_area);

%{

area_name_list= [ area_name_1 ; area_name_2 ;... ]
area_value_list = [ area_val_1 ; area_val_2]



for area_index in 1:relevant_size(area_list=[area_val_1 , area_name_1 ;area_val_2 , area_name_2 ;... ]


selected_area_1 = area_value_list(area_index, :);
%RODO: check selected_area_1 == selected_area

[s1_sa , s2_sa] = size(selected_area);
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
%imtool(mat2gray(distance_map_list_a3));



%-----------------------------------------------------------------

  
 %############## get all MSER regions




 %play with mserfeature




%MSER

%MSER parameters
min_range = 600;
max_range = 1700;
ThresholdDelta =0.4;
%,'MaxAreaVariation' ,.25);
%'ThresholdDelta',3,


%doublons filter para
%MserRegions = regions_M_242_a2;
%MserRegions = MergeM_a1;
[bg_row,bg_col,~] = size(background);
%max_aera = maxAereaRange;
slide_by_raw = bg_row; 
slide_by_col= bg_col;
min_dist=20;
plot_f_1 = false;
bg=background;
compress_percentage = .25;
same_size_percentage = .85; %.6?
min_range = min_range;
max_range = max_range;
close_enough = 20; % a reduire en tant que premier filtre qui enleve les pures doublons. genre autour de 5-10



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


concf_MSER_after_inter = MSERRegions([r1f.PixelList;r2f.PixelList;r3f.PixelList; rHf.PixelList; rEf.PixelList]);
%concatenated_MSER_regions = MSERRegions([MserRegions_1.PixelList;MserRegions_2.PixelList;MserRegions_3.PixelList; MserRegions_H.PixelList; MserRegions_E.PixelList]);



%PLOTTING

%all_3_regr
%{
figure('Name','a3_mser_regions_first_cluster');imshow(background); hold on;
hp = impixelinfo;

%figure('Name','MSER on RGB','NumberTitle','off');imshow(rgb); hold on;
%plot(regions);
%plot(regions(1),'showPixelList', true);
%plot(regions,'showEllipses',true);
plot(a1_mser_regions_first_cluster,'showPixelList',false,'showEllipses',true);


%plot(a2_mser_regions_first_cluster,'showPixelList',false,'showEllipses',true);
%plot(a3_mser_regions_first_cluster,'showPixelList',false,'showEllipses',true);


%}


%-----------------------------------------------------------------


%##############  Filter_1 : Doublons: removing closes/included ellipses

[~,regions_after_inter_filter] = filter_doublon_inter(concf_MSER_after_inter ,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );
%~~~~~~~~ ie F

%-----------------------------------------------------------------


%##############   old schol MEASURE before cell_Segm_fine_tuning_para fctn.
%{
%HnE consistency of the regions
%HnE filter check

%H_E_quantile = .8;
Area_under_PR_curve_list_by_iou_th = [];
precision_by_iou_th_list = [];
recall_by_iou_th_list = [];
for iou_treshold= 0:.05:1
    fprintf('iou th %d for area %s in tissus %s\n', iou_treshold, selected_area_name, tissue_name);

    precision_by_hne_treshold_list =[];
    recall_by_hne_treshold_list =[];
    for HnE_quantile_treshold=0:.05:1

        regions_candidates_after_hne_filter_pixel_list = candidates_mser_regions_after_hne_filter(full_H, full_E,regions_after_inter_filter, HnE_quantile_treshold, selected_area);

        [s1,s2]= size(regions_candidates_after_hne_filter_pixel_list); 
        if s1 ~= 0
             regions_candidates_after_hne_filter = MSERRegions(regions_candidates_after_hne_filter_pixel_list);      
        end



        %Plot
        %{

        figure('Name','all 3 before any filter');imshow(bg); hold on;
        hp = impixelinfo;
        plot(regions,'showPixelList',false,'showEllipses',true);
        hold off;



        figure('Name','filter duplicate of isolated layer  .......');imshow(bg); hold on;
        hp = impixelinfo;
        plot(f1_mserR_cluster1_a1,'showPixelList',false,'showEllipses',true);
        plot(f1_mserR_cluster1_a2,'showPixelList',false,'showEllipses',true);
        plot(f1_mserR_cluster1_a3,'showPixelList',false,'showEllipses',true);
        plot(Hregions,'showPixelList',false,'showEllipses',true);
        plot(Eregions,'showPixelList',false,'showEllipses',true);



        figure('Name','after all 3 duplicate filtered');imshow(bg); hold on;
        hp = impixelinfo;
        plot(candidates_regions_after_duplicate_filter,'showPixelList',false,'showEllipses',true);


        figure('Name','after hne filter');imshow(bg); hold on;
        hp = impixelinfo;
        plot(relevant_candidates_regions,'showPixelList',false,'showEllipses',true);



        %}

        %[~, rregions] =  filter_1_mser_v5(regions_candidates_after_hne_filter,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,bg_row,bg_col,slide_by_raw,slide_by_col,min_range,max_range,close_enough );



        gt_list = get_gt_list_from_folder_path (folder_path);
        pred_list = get_pred_list_from_mser_regions(regions_candidates_after_hne_filter);




        %iou_treshold= .3;
        %TODO : generique pr englober diff area et diff tissue...

        [tp , fp, fn] = get_confusion_matrix_values(pred_list,gt_list , iou_treshold);
        %ici updater+= les tp fp fn des autre area/tissue name etc labelisé
        precision = tp / (tp+fp);   %precision = tp / P
        recall = tp/(tp+fn);        %recall = tp / T

        precision_by_hne_treshold_list = [precision_by_hne_treshold_list , precision];
        recall_by_hne_treshold_list =[recall_by_hne_treshold_list, recall];

    end
    %plot( recall_by_hne_treshold_list, precision_by_hne_treshold_list);
    precision_by_iou_th_list = [precision_by_iou_th_list ;  precision_by_hne_treshold_list];
    recall_by_iou_th_list = [recall_by_iou_th_list ; recall_by_hne_treshold_list];
    Area_under_PR_curve = trapz(precision_by_hne_treshold_list,recall_by_hne_treshold_list);
    Area_under_PR_curve_list_by_iou_th = [Area_under_PR_curve_list_by_iou_th ;  Area_under_PR_curve];
end

%plot( recall_by_hne_treshold_list, precision_by_hne_treshold_list);
%dici on peut deja fixer what's the best iou_treshold
plot(recall_by_iou_th_list(1,:) , precision_by_iou_th_list(1,:))
title('Precision-Recall graph for different IoU threshold')
hold on;
[ss,~] = size(precision_by_iou_th_list);
for i=2:ss
    plot(recall_by_iou_th_list(i,:) , precision_by_iou_th_list(i,:))
end
hold off;




% 11 pts interp. for mAP
% AUC  for more para tunings.







%show_crossing(bg, pred_list, gt_list);

q=5;
%}

%############# vizu label & predict under same bg area.
%show_crossing(bg, pred_list, gt_list);

%##############  LABELIZATION.
%{
%manual labelization : 
        %struct_path = '/Users/lior/Desktop/Image & Analysis /' +  string(tissue_name) + '.mat';
        %selected_area  = [500 1000  (1000-500) (1500-1000)];
        %selected_area_name = 'area_1';
        
        
        tissue_name = 'C161VA_x20_X700Y600';
        struct_path = 'C161VA_x20_X700Y600_smoothed_and_wl_shaped.mat';



        data = load('C161VA_x20_X700Y600_smoothed_and_wl_shaped.mat');
        lambda_mat = data.lambda;
        wl_mat = data.spec;
        rgb = data.RGB;
        %imtool(rgb);
        
        [full_H, full_E] = get_HnE(wl_mat,lambda_mat);


        crop_area_1 = [1700 700 300 300];
        crop_area_2 = [1200 1500 300 200];
        crop_area_3 = [500 800 150 300];
        crop_area_4 = [800 550 400 150];



        selected_area = crop_area_4;
        selected_area_name = 'area_4';
        img = imcrop(rgb, selected_area);

        H = imcrop(full_H, selected_area);
        E = imcrop(full_E, selected_area);

        %manual_cell_labelization(tissue_name, crop_area, crop_area_name);
        manual_cell_labelization(tissue_name, img, selected_area, selected_area_name, H, E);
%}
        


%-----------------------------------------------------------------


%brouillon & old filters
%{



%[~, f1_mserR_cluster1_a1] =  f1m6(MserRegions_1,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,a,b,slide_by_raw,slide_by_col,min_range,max_range,close_enough );
%min_dist_ontry2 = filter_1_mser_v3(try2,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,a,b,max_area_raw,max_area_col,min_range,max_range );


%m_r = get_matching_rate(cell_pred_1, cell_gt_1);

%plot(mser_regions_M_242_a2,'showPixelList',false,'showEllipses',true);


[l,~] = size(f1_mserR_cluster1_a1.PixelList);
pred_all = liste sans duplicate.
for i = 1:l
    pred_all = pred_all .add( f1_mserR_cluster1_a1.PixelList(i));
end

%gt_cell_1_binary_mat =load(struct_path + area_chosen + '/cell_1.mat');
%?? gt_cell_1_binary_mat =load('/Users/lior/Desktop/Image & Analysis /cells_label/'+ string(tissue_name)+'/area_4/cell_1.mat');
for j=1:len(cells_in_selected_area) (chez moi =5 ou 2 i.e. list de tagged cell in specific area) 
gt_all_cell_union = gt_all_cell_union | gt_cell_j_binary_mat;
%cell_segm_get_score_(struct_path,pred_all, gt_all_cell_union)
%}



%{
TODO: dans f1, lorsque 2 mser proche, choisir celui qui a la plus haute moyenne de max(h_hit, E_hit)  (ccl : ca va donner du max(max(…))   ) au lieu de choisir juste le plus grand des 2.

Check why le petit en bas a droite nest pas pris (peut etre thx to 1 ca v etre ok )


%}

%try to vizualise pred and g on same bg together. need to convert gt first.



% tissues: x40_r1_r2,  P055A_X2000Y0_R1 / R2  , C083hA , C161VA_x20_X700Y600_20
            % 3.rawstack.mat , C083HA_0_1800_0
            
%}
