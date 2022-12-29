 %{
tissue_name ='x40_r1_r2';

seleted_area_list = [
 
    crop_area_1 = [500 1000  (1000-500) (1500-1000)];
    crop_area_2 = [1750 580  (2150-1750) (820-580)];
    crop_area_3 = [ 1500 1000  (2000-1500) (1500-1000)];
    crop_area_4 = [ 3400 1000  (3800-3400) (1400-1000)];
    all = [];
                    ]
    selected_area =crop_area_4 =seleted_area_list(4)
    selected_area_name = 'area_4';
    %selected_area = all;


    relevant_cluster = [0.21, 0.05, 0.37];   here for tissue: tissue_name ='x40_r1_r2'
    weighted_kmean_name = 'x24wkm_20_kmeans_and_D.mat'

%}
function [data , lambda_mat , wl_mat ,rgb,  distance_map_list] = initialization_for_one_tissue(tissue_name,weighted_kmean_name, bottom_wl_treshold, top_wl_treshold, number_of_clusters, relevant_clusters); 


    %-----------------------------------------------------------------


    %############## I initialization -uncoment

    % load tissue we work on 

    %wetransfer_3-rawstack-jpg_2021-08-30_1417/
%    tissue_name ='x40_r1_r2';
    %TO CHECK NEXT LINE 
    struct_path = pwd  + "/" + string(tissue_name) + ".mat";
    %struct_path = '/Users/lior/Desktop/Image & Analysis /' +  string(tissue_name) + '.mat';
    data = load(struct_path);
    %imtool(data.rgb)
    lambda_mat = data.lambda;
    rgb = data.rgb;
    wl_mat = data.spec;
    [a,b, ~ ] = size(wl_mat);
    
    %-----------------------------------------------------------------

    
    [~,s_name_wkm]= size(weighted_kmean_name);
    if   s_name_wkm>1
        %load weighted kmeans 
        wkm= load(weighted_kmean_name);
        wd = wkm.D;
    else
            %##############  Cutting and smothing  wl from spec

        %TOdO : filter out of wl_cut + generic filter instead of medfilt2

        %  cut and same smoothing to all WL , and  save ?  the new spectrum

        smoothed_and_wl_shaped = medfilt2(read_unique_wl(wl_mat , lambda_mat, bottom_wl_treshold));
        for i = bottom_wl_treshold+1 :top_wl_treshold
        single_wl = read_unique_wl(wl_mat , lambda_mat, i);
        smoothed_wl= medfilt2(single_wl);
        smoothed_and_wl_shaped = cat(3,smoothed_and_wl_shaped, smoothed_wl);
        end
        lambda_mat = lambda_mat(bottom_wl_treshold:top_wl_treshold);


        %SAVING: 
        
        %smoothed_and_wl_shaped  = cast(new_mat, 'double');
        %name = x40_r1_r2_spec+ '_smoothed_and_wl_shaped';
        %save(name,'smoothed_and_wl_shaped');

        %vizualization
        %{
        cutted_and_smoothed_spectral_img =  spectoRGB(smoothed_and_wl_shaped,lambda_mat, 0,'name');
        cutted_and_smoothed_spectral_img_formatted =  uint8(255 * mat2gray(cutted_and_smoothed_spectral_img));
        imshow(cutted_and_smoothed_spectral_img_formatted);
        %}
        %##############  II kmeans

         % compute weighted kmeans
        input_spectral_img = smoothed_and_wl_shaped;
        %input_spectral_img = imcrop(input_spectral_img,[2761 686  (3781-2761) (1409-686)]);

        save_kmean_img=0; %1 if want to save img
        name_of_input_img = tissue_name +'_wkm';
     %TODO: full path as new last param in wkm fct
        [W_kmean_mat_output, wd]= weighted_kmeans_algo(input_spectral_img, number_of_clusters, save_kmean_img, name_of_input_img,pwd);
        
        %{
        %imtool(wkm.kmean_mat_output);
        %km = uint8(255*(wkm.kmean_mat_output));
        %}
    end
    %-----------------------------------------------------------------
    %##############  III select good fit cluster 
    %ask user to select relevant cluster values
    % when watching %imtool(wkm.kmean_mat_output);
%    relevant_cluster = [0.21, 0.05, 0.37];
    [~ , s2]= size(relevant_clusters);
    dim_elt_in_km = [];
    for elt = 1:s2
        dim_elt_in_km = [dim_elt_in_km, round(number_of_clusters * relevant_clusters(elt))];
    end
    %----------------------------------------------------------------------
    %##############  IV DISTANCE MATRIX FROM KMEANS

    %Load, preprocess, inverse & normalize, distance of each pixel from
    %centroids and get the interesting clusters.

    %TODO: if load is true
    %load Distance 
    wd = wkm.D;
    %else : ie if not loaded and wkmean computed on line:
    %wd = D:

    %inverse and normalise:
    X =(1./wd);
    XN = normalize(X,2,'range');
    %reshape/preprocess
    [a,b]=size(wkm.kmean_mat_output);%a=2008;%b=4312;
    distance_map_list= zeros(a,b, s2);
    for elt= 1:s2  
    distance_map_list(:,:,elt) = reshape(XN(:,dim_elt_in_km(elt)), [a,b]);
    end


end