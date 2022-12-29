
% plans 
%{
-automatiser cluster selection - dimension in D wkmean
- faire schema dune une pipe line claire
- on veut du high recall avant tout . jouer avec fpr tpr pour decider best
params.
- requierd un exemple taggé

%}
%imtool(whc);


%{
%new_white =imread('C083hA_brut_Whitened.png');
whc=imcrop(new_white,[4230 115  (4960-4230) (630-115)]);
kmc=imcrop(km,[4230 115  (4960-4230) (630-115)]);
%imtool(new_white);
%}
% mask an img on another.

%{
%qa = load('aaa_255_188_215.mat');
%temp= qa.all_3_regr;
temp = all_3_regr;
mask_kernel = imcrop(temp,[4230 115  (4960-4230) (630-115)]);
origin_img = im2double(whc);
imtool(origin_img);
masked_img=add_mask_on_img(origin_img, mask_kernel);
%imtool(masked_img);
imtool(uint8(255*(masked_img)));
%}

%-----------------------------------------------------------------
%From a to z
%Load
%mat_path = '/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/C083HA_0_1800_0.mat';
%mymat = load(mat_path);
%imtool(mymat.RGB);

%{
%I whiten rgb
%{
%whitening rgb img
rgbImage = uint8(255*(mymat.RGB));
top_left_corner_x_area1 =100;
top_left_corner_y_area1= 750;
bottom_right_corner_x_area1 =560;
bottom_right_corner_y_area1 =1000;

top_left_corner_x_area2 =6400;
top_left_corner_y_area2= 0;
bottom_right_corner_x_area2 =6550;
bottom_right_corner_y_area2 =170;

save_new_white=0;
name = 'C083hA_brut';


%imtool(rgbImage);
new_white = whitening_rgb (rgbImage,top_left_corner_x_area1,top_left_corner_y_area1,bottom_right_corner_x_area1,bottom_right_corner_y_area1,top_left_corner_x_area2,top_left_corner_y_area2,bottom_right_corner_x_area2,bottom_right_corner_y_area2,save_new_white,name);
%imtool(new_white);
%}


%new_white =load('C083hA_brut_Whitened.png');
whc=imcrop(new_white,[4230 115  (4960-4230) (630-115)]);


%II kmeans
%weighted kmeans
%{
input_spectral_img = mymat.spec;
%input_spectral_img = imcrop(input_spectral_img,[2761 686  (3781-2761) (1409-686)]);
number_of_clusters = 20;
%1 if want to save img
save_kmean_img=1;
name_of_input_img = 'C083hA_brut';


[W_kmean_mat_output, D]=weighted_kmeans_algo(input_spectral_img, number_of_clusters, save_kmean_img, name_of_input_img);
imtool(W_kmean_mat_output);

%save('20_km_x40_r1_r2__smoothed_weighted_and_D','W_kmean_mat_output','D');

%}
%wkm= load('C083hA_brut_20_kmeans_and_D.mat');
%km = uint8(255*(wkm.kmean_mat_output));
%kmc=imcrop(km,[4230 115  (4960-4230) (630-115)]);

%3vals = 188 , 215, 255
% vals in D : 15 , 17 , 20


%III all 3 reg
%{
%selecting 3 clusters separatley 
%{
wkn_img=km;
% select min - max allowed size
Max_allowed_cc_size = 15000;
Min_allowed_cc_size = 20;


%select background
[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;



bg_reduction = 0;
save_filtered_img = 0;


cancerous_zone = [188];
cell_area = cancerous_zone;
cell_area_name = 'clus_188';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
f_i_188 = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
%imtool(f_i_13);


cancerous_zone = [215];
cell_area = cancerous_zone;
cell_area_name = 'clus_215';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
background = f_i_188;
bg_reduction = 1;
f_i_188_215 = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
%imtool(f_i_161_13);


cancerous_zone = [255];
cell_area = cancerous_zone;
cell_area_name = 'clus_255';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
background = f_i_188;
bg_reduction = 0;
all_3_regr = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
imtool(all_3_regr);

%save('aaa_255_188_215', 'all_3_regr');


%}
%all_3_regr = load('aaa_255_188_215.mat');
%all3_a1 =imcrop(all_3_regr,[4230 115  (4960-4230) (630-115)]);

%just 255 cluster
%{
wkn_img=km;
% select min - max allowed size
Max_allowed_cc_size = 15000;
Min_allowed_cc_size = 20;


%select background
[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;



bg_reduction = 0;
save_filtered_img = 0;



cancerous_zone = [255];
cell_area = cancerous_zone;
cell_area_name = 'clus_255';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
background = binay;
bg_reduction = 1;
filt_selec_255 = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
imtool(filt_selec_255);


%save('aaa_filt_255_only', 'filt_selec_255');
%}

%}



%IV DISTANCE MATRIX FROM KMEANS
%{
%Load, preprocess, inverse & normalize, distance of each pixel from
%centroids and get the 3 interesting clusters.

%load Distance 
wd = wkm.D;


%inverse and normalise
X =(1./wd);
XN = normalize(X,2,'range');

%reshape/ preprocess

[a,b]=size(km);
%a=2008;
%b=4312;
distance_to_Centroid_188 = reshape(XN(:,15), [a,b]);
distance_to_Centroid_215 = reshape(XN(:,17), [a,b]);
distance_to_Centroid_255 = reshape(XN(:,20), [a,b]);

distance_to_Centroid_188_a1 =imcrop(distance_to_Centroid_188,[4230 115  (4960-4230) (630-115)]);
distance_to_Centroid_215_a1 =imcrop(distance_to_Centroid_215,[4230 115  (4960-4230) (630-115)]);
distance_to_Centroid_255_a1 =imcrop(distance_to_Centroid_255,[4230 115  (4960-4230) (630-115)]);



%show output
%three_sum = crop_distance_to_Centroid_81 +crop_distance_to_Centroid_13 +crop_distance_to_Centroid_161;
%imtool(mat2gray(distance_to_Centroid_188));
%imtool(mat2gray(distance_to_Centroid_215));
%imtool(mat2gray(distance_to_Centroid_255));

%}





%V MSER
%{
min_range = 100;
max_range = 1700;
wkn_img=distance_to_Centroid_188_a1;
bg=whc;

[regions_C83_188_a1,~] = detectMSERFeatures(wkn_img, 'RegionAreaRange',[min_range ,max_range]);

wkn_img=distance_to_Centroid_215_a1;
[regions_C83_215_a1,~] = detectMSERFeatures(wkn_img, 'RegionAreaRange',[min_range ,max_range]);


wkn_img=distance_to_Centroid_255_a1;
[regions_C83_255_a1,~] = detectMSERFeatures(wkn_img, 'RegionAreaRange',[min_range ,max_range]);

%}


%VI MERGE MSER

%MERGE 3 MSER ellipses struct 
%{
%location count PixelList
PixelList = [];
for i=1:regions_C83_188_a1.Count
    PixelList =[PixelList; {regions_C83_188_a1(i).PixelList}];
end
for i=1:regions_C83_215_a1.Count
    PixelList =[PixelList; {regions_C83_215_a1(i).PixelList}];
end
for i=1:regions_C83_255_a1.Count
    PixelList =[PixelList; {regions_C83_255_a1(i).PixelList}];
end
MergeM_a1 = MSERRegions(PixelList);

%}

%PLOTTING
%{
%all_3_regr
figure('Name','merge 3');imshow(whc); hold on;
%figure('Name','MSER on RGB','NumberTitle','off');imshow(rgb); hold on;
%plot(regions);
%plot(regions(1),'showPixelList', true);
%plot(regions,'showEllipses',true);
plot(MergeM_a1,'showPixelList',true,'showEllipses',true);
%}



%VII filter: removing closes/included ellipses
%{

%filter_1_window

%option - plot the original
%{
figure();imshow(whc); hold on;
hp = impixelinfo;
plot(MergeMser);

%}

MserRegions = MergeM_a1;
%MserRegions = MergeM_a1;
[a,b,~] = size(whc);
%max_aera = maxAereaRange;
slide_by_raw = 200; 
slide_by_col= 200;
min_dist=20;
plot_f_1 = true;
bg=whc;
compress_percentage = .25;
same_size_percentage = .85; %.6?
min_range = 100;
max_range = 1700;


%new_f1_a1 =  filter_1_mser_v3(MserRegions,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,a,b,max_area_raw,max_area_col,min_range,max_range );
[short_W_PL, shortmser] =  filter_1_mser_v4(MserRegions,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,a,b,slide_by_raw,slide_by_col,min_range,max_range );
%min_dist_ontry2 = filter_1_mser_v3(try2,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,a,b,max_area_raw,max_area_col,min_range,max_range );

%}




%}

%-----------------------------------------------------------------

%general uses
%{
%km crop : kmc
%{
km = imread('x40_r1_r2_smoothed_weighted_filter_20_kmeans.png');
kmc=imcrop(km,[2761 686  (3781-2761) (1409-686)]);
%imtool(kmc);
%}
%c_km_20: == kmc from other source. 
%{
%wd = load('20_km_x40_r1_r2__smoothed_weighted_and_D.mat');
km_20 = wd.W_kmean_mat_output;
%imtool(km_20);
c_km_20 =  imcrop(km_20,[2761 686  (3781-2761) (1409-686)]);
%imtool(uint8(255*c_km_20));
%imtool(c_km_20);

%}
%rgb crop: rgbc
%{
str = load('x40_r1_r2_smoothed_weighted_filter.mat');
rgb = str.rgb;
rgbc= imcrop(rgb,[2761 686  (3781-2761) (1409-686)]);
%imtool(rgbc);
%}
%spec crop : spec_c
%{
str = load('x40_r1_r2_smoothed_weighted_filter.mat');
spec = str.smoothed_x40;
ymin = 2761;
ymax = 3781;
xmin = 686;
xmax=1409;
spec= spec(xmin+1:xmax-1,ymin+1:ymax-1,:);
%spec_c = imcrop(spec,[2761 686  (3781-2761) (1409-686)]);
%}
%wh crop : whc
%{
str = load('x40_r1_r2_smoothed_weighted_filter.mat');
wh = str.white_rgb;
whc = imcrop(wh,[2761 686  (3781-2761) (1409-686)]);
%}


%all_3_regr 
%{

%selecting 3 clusters separatley 
%wkn_img=wkn_img; 
wkn_img=kmc;
% select min - max allowed size
Max_allowed_cc_size = 150000;
Min_allowed_cc_size = 20;


%select background
[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;



bg_reduction = 0;
save_filtered_img = 0;


cancerous_zone = [13];
cell_area = cancerous_zone;
cell_area_name = 'clus_13';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
f_i_13 = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
%imtool(f_i_13);



cancerous_zone = [161];
cell_area = cancerous_zone;
cell_area_name = 'clus_161';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
background = f_i_13;
bg_reduction = 1;
f_i_161_13 = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
%imtool(f_i_161_13);


cancerous_zone = [81];
cell_area = cancerous_zone;
cell_area_name = 'clus_81';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
background = f_i_161_13;
bg_reduction = 1;
all_3_regr = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
imtool(all_3_regr);
%}
%wk_13_81
%{
wkn_img=wkn_img; 

% select min - max allowed size
Max_allowed_cc_size = 150000;
Min_allowed_cc_size = 20;


%select background
[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;



bg_reduction = 0;
save_filtered_img = 0;


cancerous_zone = [13, 81];
cell_area = cancerous_zone;
cell_area_name = 'clus_13_81';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
wk_13_81 = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
imtool(wk_13_81);
%}
%sum 13_81_161
%{
%wk_13_81 = f_i_13 + f_i_81;
%imtool(wk_13_81);
%all_3 = wk_13_81+f_i_161;
%imtool(all_3);
%}
%Eosin Source
%{
nE = imread('current_work_e.png');
Eosin = imcrop(nE,[2761 686  (3781-2761) (1409-686)]);
%imtool(Eosin);

%filtering eosin - optional
%{
Eosin_th_120= Eosin>120;
iimg = imtool(Eosin_th_120);
set(iimg, 'Name','Eosin with TH 120');
%}

%edge - optional
%{
w = edge(q);
imtool(w);
%}

%}
%create wk_13_81 
%{
wkn_img=kmc; 
% select min - max allowed size
Max_allowed_cc_size = 150000;
Min_allowed_cc_size = 20;
%select background
[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;
%para
bg_reduction = 0;
save_filtered_img = 0;
cancerous_zone = [13, 81];
cell_area = cancerous_zone;
cell_area_name = '13_81';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
wk_13_81 = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
%imtool(wk_13_81);
%}
%}
%-----------------------------------------------------------------

%MERGE 3 MSER ellipses struct 
%{
%Old
%{
temp = load('3_list_ellipse_cancerous.mat');
f_2_regions_13= temp.el_13;
f_2_regions_81= temp.el_81;
f_2_regions_161= temp.el_161;

%location count PixelList
ww = [];
Count =  f_2_regions_13.Count(1) + f_2_regions_81.Count(1) +f_2_regions_161.Count(1);
Location=[];
Axes =[];
Orientation=[];
PixelList=[];
for i=1:f_2_regions_13.Count
    Location=[Location;  f_2_regions_13(i).Location];
    Axes =[Axes;f_2_regions_13(i).Axes ];
    Orientation=[Orientation;  f_2_regions_13(i).Orientation];
    ww =[ww; {f_2_regions_13(i).PixelList}];


end

for i=1:f_2_regions_81.Count
    Location=[Location;  f_2_regions_81(i).Location];
    Axes=[Axes;  f_2_regions_81(i).Axes];
    Orientation=[Orientation;  f_2_regions_81(i).Orientation];
    ww =[ww; {f_2_regions_81(i).PixelList}];
    

end
for i=1:f_2_regions_161.Count
    Location=[Location; f_2_regions_161(i).Location];
    Axes=[Axes;  f_2_regions_161(i).Axes];
    Orientation=[Orientation;  f_2_regions_161(i).Orientation];
    ww =[ww; {f_2_regions_161(i).PixelList}];
end

MergeMser =struct('Count',Count,'Location', Location,'Axes', Axes, 'Orientation',Orientation,'PixelList',PixelList);

MergeMser.PixelList = ww;
%}

%Load if not yet
%{
temp = load('3_list_ellipse_cancerous.mat');
f_2_regions_13= temp.el_13;
f_2_regions_81= temp.el_81;
f_2_regions_161= temp.el_161;
%}


%location count PixelList
PixelList = [];
for i=1:f_2_regions_13.Count
    PixelList =[PixelList; {f_2_regions_13(i).PixelList}];
end
for i=1:f_2_regions_81.Count
    PixelList =[PixelList; {f_2_regions_81(i).PixelList}];
end
for i=1:f_2_regions_161.Count
    PixelList =[PixelList; {f_2_regions_161(i).PixelList}];
end
MergeMser = MSERRegions(PixelList);

%}

%PLOTTING
%{
%all_3_regr
figure('Name','MSER on whc');imshow(whc); hold on;
%figure('Name','MSER on RGB','NumberTitle','off');imshow(rgb); hold on;
%plot(regions);
%plot(regions(1),'showPixelList', true);
%plot(regions,'showEllipses',true);
plot(MergeMser,'showPixelList',false,'showEllipses',true);
%}


%-----------------------------------------------------------------
%DISTANCE MATRIX FROM KMEANS

%Load, preprocess, inverse & normalize, distance of each pixel from
%centroids and get the 3 interesting clusters.
%{
%load Distance matrix
wd = load('20_km_x40_r1_r2__smoothed_weighted_and_D.mat');
wd= wd.D;

%inverse and normalise
X =(1./wd);
XN = normalize(X,2,'range');

%reshape/ preprocess

[a,b]=size(rgb);
%a=2008;
%b=4312;
distance_to_Centroid_13 = reshape(XN(:,2), [a,b]);
distance_to_Centroid_81 = reshape(XN(:,7), [a,b]);
distance_to_Centroid_161 = reshape(XN(:,13), [a,b]);

crop_distance_to_Centroid_13 =  imcrop(distance_to_Centroid_13,[2761 686  (3781-2761) (1409-686)]);
crop_distance_to_Centroid_81 =  imcrop(distance_to_Centroid_81,[2761 686  (3781-2761) (1409-686)]);
crop_distance_to_Centroid_161 =  imcrop(distance_to_Centroid_161,[2761 686  (3781-2761) (1409-686)]);


%show output
%three_sum = crop_distance_to_Centroid_81 +crop_distance_to_Centroid_13 +crop_distance_to_Centroid_161;
%imtool(mat2gray(crop_distance_to_Centroid_13));
%imtool(mat2gray(crop_distance_to_Centroid_81));
%imtool(mat2gray(crop_distance_to_Centroid_161));

%}


%-----------------------------------------------------------------
% MSER + post process

%MSER
%{
min_range = 300;
max_range = 1650;
wkn_img=all_3_regr;
bg=whc;
%imtool(wkn_img);

[regions_all_3_regr,cc] = detectMSERFeatures(wkn_img, 'RegionAreaRange',[min_range ,max_range]);

%OLD
%{
%MANUALLY SELECT region area range according to ellipse size

j=1;
restrained_regions = struct([]);
Location = [];
Axes = [];
Orientation =[];
PixelList = [];
for i=1:size(regions)
    ax = regions(i).Axes;
    taille =double(ax(1,1)) + double(ax(1,2));
    if  taille > 40 && taille <250
        fprintf('%d \n ', i);
        ax = regions(i).Axes;
        loc= regions(i).Location;
        or = regions(i).Orientation;
        Location = [Location  ; [loc(1,1),loc(1,2)]];
        Axes = [Axes ; [single(ax(1,1)),single(ax(1,2))]];
        Orientation =[ Orientation ; or];
        PixelList = [PixelList ;{regions(i).PixelList}];

        j = j+1;
    end
end

%OLD TRY Manual arearange
%{
j=1;
restrained_regions = struct([]);

for i=1:size(regions)
    ax = regions(i).Axes;
    taille =double(ax(1,1)) + double(ax(1,2));
    if  taille > 50 && taille <150
        fprintf('%d \n ', i);
        ax = regions(i).Axes;
        loc= regions(i).Location;
        or = regions(i).Orientation;
        restrained_regions(j).Location = {loc(1,1),loc(1,2)};
        %restrained_regions(j).Axes(1,1) = double(ax(1,1));
        restrained_regions(j).Axes = {single(ax(1,1)),single(ax(1,2))};
        restrained_regions(j).Orientation = or;
        restrained_regions(j).PixelList = {regions(i).PixelList};

        j = j+1;
    end
end

%save('r_r',' Location','Axes','Orientation', 'PixelList');
%}

Count = j-1;
%Convert to correct format
save('r_r','Location','Axes','Orientation', 'PixelList', 'Count');
r_r = load('r_r.mat');
re_r = MSERRegions(r_r.PixelList);

%PLOTTING
figure; imshow(all_3_regr); hold on;
%plot(re_r);
%plot(re_r,'showPixelList', true);

 %}


%PLOTTING
%{
%all_3_regr
figure('Name','MSER on RGB');imshow(bg); hold on;
%figure('Name','MSER on RGB','NumberTitle','off');imshow(rgb); hold on;
%plot(regions);
%plot(regions(1),'showPixelList', true);
%plot(regions,'showEllipses',true);
plot(regions_all_3_regr,'showPixelList',false,'showEllipses',true);
%}
%}


%POST PROCESS

%Filter_1: removing closes/included ellipses

%{
%old
%{
%OLD
%{
%filter 1: when some centroide are close , keep the biggest ellipse only
%{
minimum_distance =35;

[s_loc, ~] = size(regions.Location);
distance_matrix = zeros(s_loc,s_loc);
for i=1:s_loc
    for j=1:s_loc
        loc_i= regions(i).Location;
        x_i =double(loc_i(1,1));
        y_i =double(loc_i(1,2));
        loc_j= regions(j).Location;
        x_j =double(loc_j(1,1));
        y_j =double(loc_j(1,2));
        distance_matrix(i,j)= abs(x_i - x_j) + abs(y_i - y_j);
        
    end
end


indx_to_rm = [];
for i=1:s_loc
    for j=1:s_loc
        if i==j
            continue
        end
        if distance_matrix(i,j)< minimum_distance
            ax_i = regions(i).Axes;
            size_i =double(ax_i(1,1)) + double(ax_i(1,2));
            
            ax_j = regions(j).Axes;
            size_j =double(ax_j(1,1)) + double(ax_j(1,2));
            
            if size_i >size_j
                indx_to_rm =[indx_to_rm, j];
            else
                indx_to_rm =[indx_to_rm, i];
            end
        end
    end
end

indx_to_rm =unique(indx_to_rm);

%removing them
PixelList = [];
for i=1:size(regions)
    if  ~ismember(i, indx_to_rm)
        PixelList = [PixelList ;{regions(i).PixelList}];
    end
end

[~,nn]=size(indx_to_rm);
fprintf('%d objects removed with filter_1 /n', nn);

f_1_regions = MSERRegions(PixelList);


%PLOTTING
%{
%all_3_regr
figure('Name','filter_1 on RGB');imshow(rgb); hold on;
%figure; imshow(rgb); hold on;
plot(f_1_regions);
%plot(re_r,'showPixelList', true);
%}

%}
%generic f_1
%{
min_dist=50;
MserRegions = regions;
plot_f_1 = true;
bg=rgb;
f_1_regions = filter_1_mser(MserRegions, min_dist,plot_f_1,bg);
%}
%}


MserRegions = MergeMser;
plot_f_1 = true;
bg=rgb;
compress_percentage = .25;
same_size_percentage = .75; %.6?
min_dist =30;
f_1_regions_all = filter_1_mser_relative_version(MserRegions,plot_f_1,bg,compress_percentage,same_size_percentage,min_dist);
%}

%filter_1_window

%option - plot the original
%{
figure();imshow(whc); hold on;
hp = impixelinfo;
plot(MergeMser);

%}


MserRegions = MergeMser;
[a,b,~] = size(whc);
%max_aera = maxAereaRange;
max_area_raw = 200; 
max_area_col= 200;
min_dist=50;
plot_f_1 = true;
bg=whc;
compress_percentage = .25;
same_size_percentage = .75; %.6?


new_f1_r =  filter_1_mser_v3(MserRegions,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,a,b,max_area_raw,max_area_col );
%}



%filter_2: removing ellipses with too much void
%{
%OLD
%{
% percentage of miss to reject the element
treshold_miss = .70;
indx_to_rm = [];
nbr_elt = f_1_regions.Count;

for i=1:nbr_elt
    [pxl_grp,~] = size(f_1_regions.PixelList(i));
    %fprintf('new group %d  \t there are %d pixels in this group \n', i,pxl_grp);

    miss_pixel=0;
    for j=1:pxl_grp
        %fprintf('new pixel %d \n',j);
        pxl_of_elt = f_1_regions.PixelList(i);
        y=(pxl_of_elt(j,1));
        x=(pxl_of_elt(j,2));
        if all_3_regr(x,y) == 0
            miss_pixel = miss_pixel + 1 ;
        end
    end
    if double(miss_pixel)/pxl_grp > treshold_miss
        indx_to_rm =[indx_to_rm, i];
    end
    %size(indx_to_rm)
end

%removing them
PixelList = [];
for i=1:nbr_elt
    if  ~ismember(i, indx_to_rm)
        PixelList = [PixelList ;{f_1_regions(i).PixelList}];
    end
end
[~,nn]=size(indx_to_rm);
fprintf('%d objects removed with filter_2 /n', nn);
f_2_regions = MSERRegions(PixelList);


%PLOTTING
%{
%all_3_regr
figure('Name','filter_2 on rgb');imshow(rgb); hold on;
%figure; imshow(all_3_regr); hold on;
%plot(f_2_regions);
plot(f_2_regions,'showPixelList', true);
%}


%}


%PARAM

MserR=f_1_regions_all;
working_img = all_3_regr;
void_definition_rate = .1;% under this value, a pxl is considered as void
void_in_ellipse_acceptance_rate = .2;% keep only regions that have less than x% voids
plot_f_2= true;
showPixelList=true;
showEllipses=true;
bg =rgb;

f_2_regions_all = filter_2_mser(working_img,MserR, void_in_ellipse_acceptance_rate,plot_f_2,bg,void_definition_rate,showPixelList,showEllipses);
%}

%-----------------------------------------------------------------

%MERGE 3 MSER ellipses with all_3_regr 
%{
%{
%location count PixelList
PixelList = [];
for i=1:MergeMser.Count
    PixelList =[PixelList; {MergeMser(i).PixelList}];
end
for i=1:all_3_regr.Count
    PixelList =[PixelList; {all_3_regr(i).PixelList}];
end

MergeMser_with_3_regr = MSERRegions(PixelList);

%}

%PLOTTING
%{
%all_3_regr
figure('Name','MSER_and_3_regr on whc');imshow(whc); hold on;
%figure('Name','MSER on RGB','NumberTitle','off');imshow(rgb); hold on;
%plot(regions);
%plot(regions(1),'showPixelList', true);
%plot(regions,'showEllipses',true);
plot(MergeMser_with_3_regr,'showPixelList',false,'showEllipses',true);
%}
%}

%-----------------------------------------------------------------


% seed point extraction

%{
%SEED METHODE 0 : SHRINKAGE: bad
%{

%{
bw =bwmorph(wk_13_81, 'shrink', 'inf');
zw =bwmorph(wk_13_81, 'skel', 'inf');
%imtool(bw);
%imtool(imfill(bw,'holes'));
%}

%put output on rgb to see

new_img = rgb;
new_img(: , : , 1)  =new_img(: , : , 1)  + 255*uint8(bw);
imtool(new_img);
%}


%SEED METHODE 1
%{
SE = strel("disk",5);
seed_1 = imerode(wk_13_81,SE);
%imtool(mat2gray(seed_1));
%put output on rgb to see
new_img = rgb;
new_img(: , : , 1)  =new_img(: , : , 1)  + 255*uint8(seed_1);
hToolFig = imtool(new_img);
set(hToolFig, 'Name', 'SEED1: erode img with disk 5, on rgb');

%}


%SEED METHOD 2
%{
se = offsetstrel('ball',20,20);
q = imerode(wk_13_81,se);
imtool(mat2gray(q));
w = mat2gray(q)> 0.7;
seed_2_1 = imfill(w,'holes');

%put on rgb img to see
new_img = rgb;
new_img(: , : , 3)  =new_img(: , : , 3)  + 255*uint8(seed_2_1);
hToolFig =imtool(new_img);
set(hToolFig, 'Name', 'SEED2.1: erode img with ball 20 +TH 0.7, on rgb');

%}

%{
se = offsetstrel('ball',40,40);
q = imerode(all_3_regr,se);
imtool(mat2gray(q));
p = mat2gray(q)>.50;
seed_2_2 = imfill(p,'holes');
imtool(seed_2_2);

%put on rgb img to see
new_img = rgb;
new_img(: , : , 1)  =new_img(: , : , 1)  + 255*uint8(seed_2_2);
hToolFig =imtool(new_img);
set(hToolFig, 'Name', 'SEED2.2: erode img with ball 40 +TH 0.5, on rgb');

%}


%SEED METHOD 2 bis : 2 + RegionMaxima : imregionmaxima
%{
se = offsetstrel('ball',30,30);
SE = strel("disk",4);

q = imerode(wk_13_81,se);
q=imfill(q,'holes');
imtool(mat2gray(q));
BW = imregionalmax(q);
imtool(BW);

%}



%CC
%{
CC = bwconncomp(z);
%CC.NumObjects

%CC.PixelIdxList
stats = regionprops('Table',CC,'basic');
centroid_list = stats.Centroid;
%ctrd = centroid_list(i,:);


x =round(centroid_list(1,1));
y =round(centroid_list(1,2));
%}

%region growing
%{
%I = logical(z);
%j = regiongrowing(I,x,y,0.2 );
%imtool(j);
%}

%EXTRACTFEATURES fct
%{
[a,b] = extractFeatures(all_3_regr, regions);
%plot(b);
%cc.PixelIdxList
%}

%put on rgb img to see
%{
new_img = rgb;
new_img(: , : , 1)  =new_img(: , : , 1)  + 255*uint8(a);
imtool(new_img);
%}



%}

%-----------------------------------------------------------------


%MSER
%{
mask = seed_1;
%mask=seed_2_1;
%mask=seed_2_2;
%mask=mat2gray(Eosin_th_120);
%mask=all_3_regr;


min_range = 300;
max_range = 1650;
[regions,cc] = detectMSERFeatures(mask, 'RegionAreaRange',[min_range ,max_range]);

%PLOTTING
figure('Name','MSER on mask');imshow(rgb); hold on;
%figure('Name','MSER on RGB','NumberTitle','off');imshow(rgb); hold on;
%plot(regions);
%plot(regions(1),'showPixelList', true);
plot(regions,'showEllipses',true);
plot(f_2_regions,'showPixelList', true,'showEllipses',false);
%plot(regions,'showPixelList',true,'showEllipses',true);


%}

%GENERALIZED MSER
%{
mask = seed_1;
%mask=seed_2_1;
%mask=seed_2_2;
%mask=mat2gray(Eosin_th_120);
%mask=all_3_regr;


wkn_img = seed_2_1;
min_range = 150;
max_range = 1650;

f_1_min_distance = 25;
f_2_treshold_miss = .80;
plot_mser=false;
plot_f_1=false;
plot_f_2=true;
bg= rgb;
regions_seed_2_1 = MSER_and_post_process(wkn_img,min_range , max_range, f_1_min_distance,f_2_treshold_miss,plot_mser, plot_f_1, plot_f_2, bg);

%}



%Merging  f_2 of mask and f_2 of all_3_regr
%{
mask =seed_2_1;

%creation of both f_2

%supposed to be run
min_range = 150;
max_range = 1650;

f_1_min_distance = 25;
f_2_treshold_miss = .80;
plot_mser=false;
plot_f_1=false;
plot_f_2=true;
bg= rgb;
regions_mask = MSER_and_post_process(mask,min_range , max_range, f_1_min_distance,f_2_treshold_miss,plot_mser, plot_f_1, plot_f_2, bg);
f_2_regions= MSER_and_post_process(all_3_regr,300 , max_range, f_1_min_distance,f_2_treshold_miss,plot_mser, plot_f_1, plot_f_2, bg);


%PLOT
%{
%all_3_regr
figure('Name','seed_1 and all_3 on rgb before');imshow(bg); hold on;
%figure; imshow(all_3_regr); hold on;
plot(regions_seed_1);
plot(f_2_regions,'showPixelList', false);
%}



new_pxl_list= [regions_mask.PixelList ; f_2_regions.PixelList];

merge_regions = MSERRegions(new_pxl_list);

%generic f_1

min_dist=50;
MserRegions = merge_regions;
plot_f_1 = true;
bg=rgb;
f_1_regions = filter_1_mser(MserRegions, min_dist,plot_f_1,bg);



%filter 1: when some centroide are close , keep the biggest ellipse only
%{
minimum_distance =50;

[s_loc, ~] = size(merge_region.Location);
distance_matrix = zeros(s_loc,s_loc);
for i=1:s_loc
    for j=1:s_loc
        loc_i= merge_region(i).Location;
        x_i =double(loc_i(1,1));
        y_i =double(loc_i(1,2));
        loc_j= merge_region(j).Location;
        x_j =double(loc_j(1,1));
        y_j =double(loc_j(1,2));
        distance_matrix(i,j)= abs(x_i - x_j) + abs(y_i - y_j);
        
    end
end


indx_to_rm = [];
for i=1:s_loc
    for j=1:s_loc
        if i==j
            continue
        end
        if distance_matrix(i,j)< minimum_distance
            ax_i = merge_region(i).Axes;
            size_i =double(ax_i(1,1)) + double(ax_i(1,2));
            
            ax_j = regions(j).Axes;
            size_j =double(ax_j(1,1)) + double(ax_j(1,2));
            
            if size_i >size_j
                indx_to_rm =[indx_to_rm, j];
            else
                indx_to_rm =[indx_to_rm, i];

            end
            
        end
    end
end
indx_to_rm =unique(indx_to_rm);
size(indx_to_rm);

%removing them
PixelList = [];
for i=1:size(merge_region)
    if  ~ismember(i, indx_to_rm)
        PixelList = [PixelList ;{merge_region(i).PixelList}];
    end
end

[~,nn]=size(indx_to_rm);
fprintf('%d objects removed with filter_1 \n', nn);

f_1_regions = MSERRegions(PixelList);


%PLOTTING
%{
%all_3_regr
figure('Name','filter_1 on RGB');imshow(rgb); hold on;
%figure; imshow(rgb); hold on;
plot(f_1_regions);
%plot(re_r,'showPixelList', true);
%}

%}

%generic f_2
wkn_img= mask + all_3_regr;
treshold_miss = .80;
MserRegions=f_1_regions;
plot_f_2= true;
bg =rgb;
f_2_regions = filter_2_mser(wkn_img, MserRegions, treshold_miss,plot_f_2,bg);


%filter_2: removing ellipses that surround void
%{
% percentage of miss to reject the element
treshold_miss = .80;

indx_to_rm = [];
nbr_elt = f_1_regions.Count;
for i=1:nbr_elt
    [pxl_grp,~] = size(f_1_regions.PixelList(i));
    %fprintf('new group %d  \t there are %d pixels in this group \n', i,pxl_grp);

    miss_pixel=0;
    for j=1:pxl_grp
        %fprintf('new pixel %d \n',j);
        pxl_of_elt = f_1_regions.PixelList(i);
        y=(pxl_of_elt(j,1));
        x=(pxl_of_elt(j,2));
        if all_3_regr(x,y) == 0
            miss_pixel = miss_pixel + 1 ;
        end
    end
    if double(miss_pixel)/pxl_grp > treshold_miss
        indx_to_rm =[indx_to_rm, i];
    end
    %size(indx_to_rm)
end

%removing them
PixelList = [];
for i=1:nbr_elt
    if  ~ismember(i, indx_to_rm)
        PixelList = [PixelList ;{f_1_regions(i).PixelList}];
    end
end
[~,nn]=size(indx_to_rm);
fprintf('%d objects removed with filter_2 \n', nn);
f_2_regions = MSERRegions(PixelList);


%PLOTTING
%{
%all_3_regr
figure('Name','seed_1 and all_3 on rgb after');imshow(rgb); hold on;
%figure; imshow(all_3_regr); hold on;
plot(f_2_regions);
%plot(f_2_regions,'showPixelList', false);
%}


%}

%}




%rm MSER Region from masks individually
%{
 %TODO: work on full hit , partial etc. nned to do CC before >..?
mask = seed_1;
[a,b] = size(mask);



%CC by bwconncomp + regionprops
%{
CC = bwconncomp(mask);
%CC.NumObjects

%CC.PixelIdxList
stats = regionprops('Table',CC,'basic');
%centroid_list = stats.Centroid;
%ctrd = centroid_list(i,:);


%x =round(centroid_list(1,1));
%y =round(centroid_list(1,2));
%}

%CC to Mask by bwlabel

[L,~] = bwlabel(mask);

% count # pixels in each CC of L.
map = nbr_elt_in_clusters(a,b,L);
%{
for un= 1:a
    for deux =1:b
        if isKey(map, L(un,deux))
            mask(un,deux)=0;
        end
    end
end
 %}





% completing the function

imtool(mask);
treshold_hit = .90;
%remove object from mser from mask
PixelList = [];
full_hit = [];

nbr_obj_in_mser = f_2_regions.Count;
for i=1:nbr_obj_in_mser

    [pxl_grp,~] = size(f_2_regions.PixelList(i));
    %fprintf('new group %d  \t there are %d pixels in this group \n', i,pxl_grp);

    hit_pixel=0;
    for j=1:pxl_grp
        %fprintf('new pixel %d \n',j);
        pxl_of_elt = f_2_regions.PixelList(i);
        y=(pxl_of_elt(j,1));
        x=(pxl_of_elt(j,2));
        if  mask(x,y) ~= 0 & isKey(map, L(un,deux))
            mask(un,deux)=0;
        end
        if mask(x,y) ~= 0
            hit_pixel = hit_pixel + 1 ;
        end
        mask(x,y)= 0;
    end
    if double(hit_pixel)/pxl_grp >= treshold_hit
        full_hit =[full_hit, i];
    end

end


%f_3_regions = MSERRegions(mask);


%PLOTTING
%{
%all_3_regr
figure('Name','rm mser from mask');imshow(mask); hold on;
%figure; imshow(all_3_regr); hold on;
%plot(f_2_regions);
%plot(f_2_regions,'showPixelList', true);
%}

%}

%-----------------------------------------------------------------

%KMEANS get D 
%{
%mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/C161VA_x20_X700Y600_Whitened_2');
mymat = load("x40_r1_r2_smoothed_weighted_filter.mat");

input_spectral_img = mymat.smoothed_x40(:,:,:);
%input_spectral_img = imcrop(input_spectral_img,[2761 686  (3781-2761) (1409-686)]);
number_of_clusters = 20;
%1 if want to save img
save_kmean_img=0;
name_of_input_img = 'x_40_r1_r2';


[W_kmean_mat_output, D]=weighted_kmeans_algo(input_spectral_img, number_of_clusters, save_kmean_img, name_of_input_img);
imtool(W_kmean_mat_output);

%save('20_km_x40_r1_r2__smoothed_weighted_and_D','W_kmean_mat_output','D');

%}

%-----------------------------------------------------------------

%filtering eosin
%{
q = nE>120;
c_13_81= f_i_13 +  f_i_81;
e_13_81= edge(c_13_81);
imtool(c_13_81 + f_i_161);
%}

%cluster edge and blurr
%{
km = imread('x40_r1_r2_smoothed_weighted_filter_20_kmeans.png');
wkn_img=imcrop(km,[2761 686  (3781-2761) (1409-686)]);
imtool(wkn_img);


%selecting 3 clusters separatley 
%{
% select min - max allowed size
Max_allowed_cc_size = 150000;
Min_allowed_cc_size = 20;


%select background
%{
%kernel_out_8_300 = imread('kernel_out_8_300.png');
%support_name= 'kernel_8_300';
%support=mat2gray(kernel_out_8_300);

%m_k_8_300 = imread('kernel_8_300_on_empty_filtered_img.png');
%support_name= 'm_k_8_300';
%support=mat2gray(m_k_8_300);

%membrane_kernel_8_100 = imread('kernel_out_8_100.png');
%support =mat2gray(membrane_kernel_8_100);
%support_name= 'membrane_kernel_8_100';


[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;


background = mat2gray(imread('C161VA_x20_X700Y600_20_kmeans_k_13_20_1000_on_empty.png'));
bg_name= 'k13_20_1000';

%}
[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;



bg_reduction = 0;
save_filtered_img = 0;


cancerous_zone = [13];
cell_area = cancerous_zone;
cell_area_name = 'clus_13';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
f_i_13 = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
imtool(f_i_13);



cancerous_zone = [161];
cell_area = cancerous_zone;
cell_area_name = 'clus_161';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
f_i_161 = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
imtool(f_i_161);


cancerous_zone = [81];
cell_area = cancerous_zone;
cell_area_name = 'clus_81';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
f_i_81 = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
imtool(f_i_81);
%}

%same for hne separately
%{
img = load('x40_r1_r2_smoothed_weighted_filter.mat');
FImage = img.smoothed_x40;
lambda = img.lambda;
[Fraction,Err]=SVDAnalysisHnE(FImage,lambda,'g',0);
W1=[0,0,2,0; 0,0.5,0,0; 1,0,0,1];
W2=[0,0,3,0; 0,0.7,0,0; 1,0,0,1];

E= Fraction(:,:,2);
%imtool(E);

nE = uint8(255 * mat2gray(E));
imtool(nE);
%imwrite(nE, 'current_work_e.png');

H= Fraction(:,:,3);
%imtool(H);

nH = uint8(255 * mat2gray(H));
%imtool(nH);
%imwrite(nH, 'current_work_h.png');


Eosin = imcrop(nE,[2761 686  (3781-2761) (1409-686)]);
imtool(Eosin);
Hematoxilyn =imcrop(nH,[2761 686  (3781-2761) (1409-686)]);


function FImage=FixedPatternCorrect(Image,lambda,EmptyRef,SmoothLength,WB)

% Image - for correction
% lambda - vector of wavelengths
% EmptyRef - a reference matrix produced from an Empty screen (or part of
%       image) using FixedPatternLearn.m
% SmoothLength: e.g. =1 smoothing along Y before applying calibration with EmptyRef
% WB: 1 - perform White balance
NwlMax=size(Image,3);
MeasSkipLength=size(EmptyRef,2);       % the skip length during the measurement 
FImage=zeros(size(Image));


% Optional smoothing along Y axis - improve SNR but may hurt resolution:
if SmoothLength>0
    SmoothVector=[(1:1:SmoothLength) (SmoothLength:-1:1)];  
    SmoothVector=SmoothVector'/sum(SmoothVector);
    SmthImage=zeros(size(Image));
    for Nwl=1:NwlMax
         SmthImage(:,:,Nwl)=conv2(Image(:,:,Nwl),SmoothVector,'same');
    end
else
    SmthImage=Image;
end

% normalizing the whole image accoring to reference provided (EmptyRef):   

for Nwl=1:NwlMax     % loop over wavelengths
           
    for j=1:size(Image,2)
        k=mod(j,MeasSkipLength);        % finding the correct normalization vector
        if k==0; k=MeasSkipLength; end
        FImage(:,j,Nwl)=SmthImage(:,j,Nwl)./EmptyRef(:,k,Nwl);
    end
    
end
%{
RGBHE=Spec2RGB(FImage,lambda);

if WB==1
    FRGB1=RGBHE;
    figure
    subplot(1,2,1);
    imshow(FRGB1)
    h=gca;
    h.Title.String='mark white area';
    rect=round(getrect);

    for c=1:3
        WhiteBal(c)=mean(mean(FRGB1(rect(2)+[1:rect(4)],rect(1)+[1:rect(3)],c)));
        RGBHE(:,:,c)=FRGB1(:,:,c)*(1/WhiteBal(c));
    end
    
    subplot(1,2,2);
    imshow(RGBHE)
end
%}

end
function EmptyRef=FixedPatternLearn(RefImage,MeasSkipLength,YSmoothLength)

% MeasSkipLength=14;       % skip length along X during the measurement 
% YSmoothLength=12;        % half length for smoothing along Y

NwlMax=size(RefImage,3);

SmoothVector=[(1:1:YSmoothLength) (YSmoothLength:-1:1)];  
SmoothVector=SmoothVector'/sum(SmoothVector);

Simage=zeros(size(RefImage));
EmptyRef=zeros(size(RefImage,1),MeasSkipLength,NwlMax);

% finding normalization matrix:

for Nwl=1:NwlMax     % loop over wavelengths
    
    % average along Y axis - improve SNR if no info in the figure (e.g. empty slide or empty part of figure)
    Simage(:,:,Nwl)=conv2(RefImage(:,:,Nwl),SmoothVector,'same');
    
    % average in jumps equivalent to X skipping.  notice this is done for full
    % vectors hence for each raw along the X scanning direction data is
    % collected separately in order to deal with nonuniformity along y axis.
    for StartCol=1:MeasSkipLength
        EmptyRef(:,StartCol,Nwl)=mean(Simage(:,(StartCol:MeasSkipLength:size(RefImage,2)),Nwl),2);
    end
       
end
end
function [Fraction,Err]=SVDAnalysisHnE(FImage,lambda,ProgRep,Ssmooth)

% ProgRep='t';         % type of progress report:  t=text, g=graph, other= none
% Ssmooth=1;           % Smoothing on WL basis befor analysis, the values is the sigma of a gaussian convoluted with image. can take any value. 0= no smoothing. 

%===================================
% SVD analysis
%===================================
load('REF_Absorption_Spectra.mat')
clear REF_Spec LA 

Absorption=-log10(FImage);

% PARAMETERS:
% focus data on interesting part ofg spectrum, e.g. (500-700nm):
LambdaMin=475;      % min/max wavelengths of interest for the analysis
LambdaMax=700;

Ns=4;                % number of input spectra, including Hematoxylin, Eosin and two for background
dn=1;                % sub-sampling in X,Y  to reduce amount of computation by dn^2, if needed
Nr=50000;            % progress report rate
%...........................................................................

Ny=size(Absorption,1);
Nx=size(Absorption,2);
NL=size(Absorption,3);


LIdx=(lambda>LambdaMin).*(lambda<LambdaMax);
Lmin=find(LIdx, 1 ); 
Lmax=find(LIdx, 1, 'last' );             
L=lambda(Lmin:Lmax)';
NL=Lmax-Lmin+1;          


REF_Spec(:,1)=L;  %wavelength
% Actuall absorption curves of Eosin:
REF_Spec(:,2)=interp1(Eosin_Abs_Rel(:,1),Eosin_Abs_Rel(:,2),REF_Spec(:,1));
% Gaussian approximaltion to Eosin:
% WLm2=520;  Sig2=10; 
% REF_Spec(:,2)=exp(-(REF_Spec(:,1)-WLm2).^2/(2*Sig2^2));

% Actuall absorption curves of Hematoxylin:
REF_Spec(:,3)=interp1(Hematoxylin_Abs_Rel(:,1),Hematoxylin_Abs_Rel(:,2),REF_Spec(:,1));
% Gaussian approximaltion to Hematoxylin:
% WLm3=590;  Sig3=40; 
% REF_Spec(:,3)=exp(-(REF_Spec(:,1)-WLm3).^2/(2*Sig3^2));

% Additional curves to account for background and unknown materials
WLm4=700;  Sig4=50; 
REF_Spec(:,4)=exp(-(REF_Spec(:,1)-WLm4).^2/(2*Sig4^2));
WLm5=460;  Sig5=30; 
REF_Spec(:,5)=exp(-(REF_Spec(:,1)-WLm5).^2/(2*Sig5^2));

% figure(5); plot(REF_Spec(:,1),REF_Spec(:,2),REF_Spec(:,1),REF_Spec(:,3))
% .............................................

A=REF_Spec(:,2:Ns+1);               % This is the transfer matrix from fractions to combined spectrum
[U,S,V] = svd(A);                        % SVD decomposition
S1=[inv(S(1:Ns,1:Ns)) zeros(Ns,NL-Ns)];                % this is the "inverse" of S

UnMix=inv(V')*S1*inv(U);            % this is generalized inv(A), allowing to go from spectra to fractions

%========================================================================================
% PREFILTERING:

LA=zeros(Ny,Nx,NL);

% Optional smoothing of the data (each WL seperately) in (X,Y) before analysis :
for j=1:NL
    if Ssmooth>0
        LA(:,:,j)=SpatialFilter(Absorption(:,:,Lmin+j-1),Ssmooth);    
    else
        LA(:,:,j)=Absorption(:,:,Lmin+j-1);
    end
end

% Optional subsampling:
if dn~=1   
    LA=LA(1:dn:end,1:dn:end,:);
end

%========================================================================================
% Running SVD analysis per spectrum:

Fraction=zeros(Ny,Nx,Ns);
Err=zeros(Ny,Nx);

n=0;
Tstart=now;

for j=1:Ny
    for k=1:Nx
        n=n+1;
        S=squeeze(LA(j,k,:));
        fr=UnMix*S;
        Fraction(j,k,:)=fr;                    % Fraction of each component is the output of analysis
        Err(j,k)=std(S-A*fr)^2;
        
        if fix(n/Nr)==n/Nr                     % progress report
            pFinish=100*n/(Nx*Ny);
            tPass= (now-Tstart)*24*60;
            tLeft = tPass*(100-pFinish)/pFinish;
            if ProgRep=='t'
                % progress report by text:
                disp([num2str(pFinish) ' % ||  Time from start: ' num2str(tPass ) ' min. ||  Remaining Time: ' num2str(tLeft)  ' min.  '])
            elseif ProgRep=='g'
                % progress report by graph:
                figure(22); clf;
                C1=A(:,1)*fr(1);
                C2=A(:,2)*fr(2);
                C3=A(:,3)*fr(3);
                C4=A(:,4)*fr(4);
                plot(L,S,'or', L,A*fr,'-b',L,C1,'--c',L,C2,'--c',L,C3,'--c',L,C4,'--c')
                drawnow
                title(['Error (SSE): ' num2str( Err(Ny,Nx)) ' ||  Finished: ' num2str(pFinish) ' % ||  Time from start: ' num2str(tPass ) ' min. ||  Remaining Time: ' num2str(tLeft)  ' min.  ']);
                pause(0.1)
            end
        end
    end
end
end
function [IndAll,Spectra]=ShowSpec2(RGBHE,spec,lambda,n,Fraction)
% RGBHE - RGB Image for the original H&E transmission
% spec - data cube
% labmbda ...
% n - number of nearest neibors for averaging when displaying spectrum.

FRGB=FalseFractionRGB(Fraction,[0,0,3,0; 0,0.7,0,0; 1,0,0,1]);      % false RGB representing SVD analysis results (Fraction)

f=figure(6);
if isempty(f.CurrentAxes)               % if figure exist no need to redraw so magnif remains as before
    subplot(2,2,1); 
    imshow(RGBHE)  
    a=f.CurrentAxes;
    subplot(2,2,2);
    imshow(FRGB)
    b=f.CurrentAxes;
end

a.Position=[0.0    0.48    0.6    0.5];
b.Position=[0.5    0.48    0.6    0.5];

IndAll=[];
Spectra=[];

stop=0;

while stop~=1
    [x,y,button]=ginput(1);
    X=round(x);
    Y=round(y);
    Ind=[X Y];
    if button=='s' || button=='S'       % 's' for stopping
        stop=1; 
    elseif button=='p' || button=='P'   % 'p' may be used for pausing and changing magnification, hit Enter when done 
        pause
        % adjust magnif similar for both views:
        c=f.CurrentAxes;    % read imag limits from last touched subplot
        XL=c.XLim;          
        YL=c.YLim;
        a.XLim=XL;          % apply same limits to both plots
        a.YLim=YL;
        b.XLim=XL;
        b.YLim=YL;
    elseif button=='n' || button=='N'   % 'n' for changing the averaging area
        n=input('Averaging area size? : ');
    else                                % mouse-click - choosing a point a drawing spectrum
        subplot(2,2,3); 
        a3=f.CurrentAxes;
        a3.Position=[0.1300    0.0800    0.3347    0.3];
        if n==0
            Spectrum=squeeze(spec(Y,X,:));
            F1=Fraction(Y,X,1);
            FE=Fraction(Y,X,2);
            FH=Fraction(Y,X,3);
            F4=Fraction(Y,X,4);
        else
            Spectrum=squeeze(mean(spec(Y-n:Y+n,X-n:X+n,:),[1 2]));
            F1=mean(Fraction(Y-n:Y+n,X-n:X+n,1),[1 2]);
            FE=mean(Fraction(Y-n:Y+n,X-n:X+n,2),[1 2]);
            FH=mean(Fraction(Y-n:Y+n,X-n:X+n,3),[1 2]);           
            F4=mean(Fraction(Y-n:Y+n,X-n:X+n,4),[1 2]);           
        end
        plot(lambda,Spectrum)
        axis([400 800 0 1.2])
        title(['X:  ' num2str(X) ' , Y:  ' num2str(Y)])
        IndAll=[IndAll; Ind];
        Spectra=[Spectra; Spectrum];
        
        subplot(2,2,4)
        a4=f.CurrentAxes;
        a4.Position=[0.5703    0.1100    0.3347    0.3412];

        bar(1:4,[F1 FE FH F4])
        axis([0 5 -.2 1.2])
        title(['F1:  ' num2str(F1,'%.3f')   '      , Eosin:  ' num2str(FE,'%.3f') '       ,  Hematoxylin:  ' num2str(FH,'%.3f') '      , F4:  ' num2str(F4,'%.3f') ])
        text(1,1,['n= ' num2str(n)])
  
    end
   
end

end
function Af=SpatialFilter(A,Nf)

% Nf - filter half-size:
x=-Nf:Nf;
y=x';
X=ones(2*Nf+1,1)*x;
Y=y*ones(1,2*Nf+1);
F=exp(-2*(X.^2+Y.^2)/Nf^2);        % Gaussian filter definition
F=F/sum(sum(F));         

Af=conv2(A,F,'same');        % applying filter to Mask
end
function FRGB=FalseFractionRGB(Fraction,W)

% W = weights for RGB false coloring
% e.g. W=[0,0,2,0; 0,0.5,0,0; 1,0,0,1]


FRGB(:,:,1)=W(1,1)*Fraction(:,:,1)+W(1,2)*Fraction(:,:,2)+W(1,3)*Fraction(:,:,3)+W(1,4)*Fraction(:,:,4);      % R
FRGB(:,:,2)=W(2,1)*Fraction(:,:,1)+W(2,2)*Fraction(:,:,2)+W(2,3)*Fraction(:,:,3)+W(2,4)*Fraction(:,:,4);      % G
FRGB(:,:,3)=W(3,1)*Fraction(:,:,1)+W(3,2)*Fraction(:,:,2)+W(3,3)*Fraction(:,:,3)+W(3,4)*Fraction(:,:,4);      % B
end
function XRGB=Spec2RGB(spec,lambda)

R=find(lambda>600);
G=find((lambda>500).* (lambda<600));
B=find(lambda<500);
Dim=size(spec);

XRGB(:,:,1)=sum(spec(:,:,R),3);
Map=squeeze(XRGB(:,:,1));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,1)=XRGB(:,:,1)./prctile(Sig,98);

XRGB(:,:,2)=sum(spec(:,:,G),3);
Map=squeeze(XRGB(:,:,2));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,2)=XRGB(:,:,2)./prctile(Sig,98);

XRGB(:,:,3)=sum(spec(:,:,B),3);
Map=squeeze(XRGB(:,:,3));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,3)=XRGB(:,:,3)./prctile(Sig,98);
end
%}

%blurring
%{
f_i = f_i_13 + f_i_81;
o = ones(10,10)/(10*10);
b = conv2(double(f_i),o, 'same');
imtool(mat2gray(b/256));
%}

%Threshold the Image
%{
bw = imbinarize(b);

%fill gap 
se = strel('disk',2);
bw = imclose(bw,se); 
bw = imfill(bw,'holes');
imshow(bw);

%}


%edg  = edge(b/256, 'Sobel');
%imtool(edg);

%}

%edge and blurr old
%{
%{
rgb = imread("x40_r1_r2_smoothed_Whitened.png");

%imtool(rgb);
rgb=imcrop(rgb,[2761 749  (3781-2761) (1409-749)]);



e = imread('current_work_e.png');
e=imcrop(e,[2761 749  (3781-2761) (1409-749)]);
%imtool(e);
q = e>120;
e_edge = edge(q,'Sobel');
%imtool(e_edge);


h= imread('current_work_h.png');
h=imcrop(h,[2761 749  (3781-2761) (1409-749)]);
%imtool(h);

c_255 = imread('current_work_clus_255.png');
c255_edge = edge(c_255,'Sobel');
%imtool(c255_edge);

c_161 = imread('current_work_clus_161.png');
c161_edge = edge(c_161,'Sobel');
%imtool(c161_edge);
c_54 = imread('current_work_clus_54.png');
c54_edge = edge(c_54,'Sobel');
%imtool(c54_edge);

all_C = uint8(c_161) + uint8(c_255) + uint8(c54_edge);
%imtool(all_C);
%}

mask_kernel= uint8(all_C);
origin_img = rgb;
masked_img=add_mask_on_img(origin_img, mask_kernel);
%imtool(masked_img);


o = ones(15,15)/(15*15);

b = conv2(double(all_C),o, 'same');
%imtool(b/256);


edg  = edge(b/256, 'Sobel');
imtool(edg);
%}
%-----------------------------------------------------------------

%current WORK: 3 img saved of smoothed x40 - last wl , 20kmeans , clusters of
%cancerous separately: 54 161 255
%filtered_selection ? 
%{
%filtered_selection

%parameters

wkn_img = imread('x40_r1_r2_smoothed_meanfilter_5_20_kmeans.png');
wkn_img=imcrop(wkn_img,[2761 749  (3781-2761) (1409-749)]);

% select min - max allowed size
Max_allowed_cc_size = 15000000;
Min_allowed_cc_size = 20;


%select background
%{
%kernel_out_8_300 = imread('kernel_out_8_300.png');
%support_name= 'kernel_8_300';
%support=mat2gray(kernel_out_8_300);

%m_k_8_300 = imread('kernel_8_300_on_empty_filtered_img.png');
%support_name= 'm_k_8_300';
%support=mat2gray(m_k_8_300);

%membrane_kernel_8_100 = imread('kernel_out_8_100.png');
%support =mat2gray(membrane_kernel_8_100);
%support_name= 'membrane_kernel_8_100';


[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;


background = mat2gray(imread('C161VA_x20_X700Y600_20_kmeans_k_13_20_1000_on_empty.png'));
bg_name= 'k13_20_1000';

%}
[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;



bg_reduction = 0;
save_filtered_img = 0;


cancerous_zone = [54];
cell_area = cancerous_zone;
cell_area_name = 'clus_54';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
f_i = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
imtool(f_i);



cancerous_zone = [161];
cell_area = cancerous_zone;
cell_area_name = 'clus_161';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
f_i = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
imtool(f_i);


cancerous_zone = [255];
cell_area = cancerous_zone;
cell_area_name = 'clus_255';
name = 'current_work';
final_name = append(name, '_', cell_area_name);
f_i = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,8);
imtool(f_i);
%}

%same for hne separately + hne functions
%{
img = load('x40_r1_r2_smoothed_meanfilter_5.mat');
FImage = img.spec;
lambda = img.lambda;
[Fraction,Err]=SVDAnalysisHnE(FImage,lambda,'g',0);
W1=[0,0,2,0; 0,0.5,0,0; 1,0,0,1];
W2=[0,0,3,0; 0,0.7,0,0; 1,0,0,1];

E= Fraction(:,:,2);
%imtool(E);

nE = uint8(255 * mat2gray(FRGB));
%imtool(nE);
%imwrite(nE, 'current_work_e.png');

H= Fraction(:,:,3);
%imtool(H);

nH = uint8(255 * mat2gray(FRGB));
imtool(nH);
%imwrite(nH, 'current_work_h.png');

function FImage=FixedPatternCorrect(Image,lambda,EmptyRef,SmoothLength,WB)

% Image - for correction
% lambda - vector of wavelengths
% EmptyRef - a reference matrix produced from an Empty screen (or part of
%       image) using FixedPatternLearn.m
% SmoothLength: e.g. =1 smoothing along Y before applying calibration with EmptyRef
% WB: 1 - perform White balance
NwlMax=size(Image,3);
MeasSkipLength=size(EmptyRef,2);       % the skip length during the measurement 
FImage=zeros(size(Image));


% Optional smoothing along Y axis - improve SNR but may hurt resolution:
if SmoothLength>0
    SmoothVector=[(1:1:SmoothLength) (SmoothLength:-1:1)];  
    SmoothVector=SmoothVector'/sum(SmoothVector);
    SmthImage=zeros(size(Image));
    for Nwl=1:NwlMax
         SmthImage(:,:,Nwl)=conv2(Image(:,:,Nwl),SmoothVector,'same');
    end
else
    SmthImage=Image;
end

% normalizing the whole image accoring to reference provided (EmptyRef):   

for Nwl=1:NwlMax     % loop over wavelengths
           
    for j=1:size(Image,2)
        k=mod(j,MeasSkipLength);        % finding the correct normalization vector
        if k==0; k=MeasSkipLength; end
        FImage(:,j,Nwl)=SmthImage(:,j,Nwl)./EmptyRef(:,k,Nwl);
    end
    
end
%{
RGBHE=Spec2RGB(FImage,lambda);

if WB==1
    FRGB1=RGBHE;
    figure
    subplot(1,2,1);
    imshow(FRGB1)
    h=gca;
    h.Title.String='mark white area';
    rect=round(getrect);

    for c=1:3
        WhiteBal(c)=mean(mean(FRGB1(rect(2)+[1:rect(4)],rect(1)+[1:rect(3)],c)));
        RGBHE(:,:,c)=FRGB1(:,:,c)*(1/WhiteBal(c));
    end
    
    subplot(1,2,2);
    imshow(RGBHE)
end
%}

end
function EmptyRef=FixedPatternLearn(RefImage,MeasSkipLength,YSmoothLength)

% MeasSkipLength=14;       % skip length along X during the measurement 
% YSmoothLength=12;        % half length for smoothing along Y

NwlMax=size(RefImage,3);

SmoothVector=[(1:1:YSmoothLength) (YSmoothLength:-1:1)];  
SmoothVector=SmoothVector'/sum(SmoothVector);

Simage=zeros(size(RefImage));
EmptyRef=zeros(size(RefImage,1),MeasSkipLength,NwlMax);

% finding normalization matrix:

for Nwl=1:NwlMax     % loop over wavelengths
    
    % average along Y axis - improve SNR if no info in the figure (e.g. empty slide or empty part of figure)
    Simage(:,:,Nwl)=conv2(RefImage(:,:,Nwl),SmoothVector,'same');
    
    % average in jumps equivalent to X skipping.  notice this is done for full
    % vectors hence for each raw along the X scanning direction data is
    % collected separately in order to deal with nonuniformity along y axis.
    for StartCol=1:MeasSkipLength
        EmptyRef(:,StartCol,Nwl)=mean(Simage(:,(StartCol:MeasSkipLength:size(RefImage,2)),Nwl),2);
    end
       
end
end
function [Fraction,Err]=SVDAnalysisHnE(FImage,lambda,ProgRep,Ssmooth)

% ProgRep='t';         % type of progress report:  t=text, g=graph, other= none
% Ssmooth=1;           % Smoothing on WL basis befor analysis, the values is the sigma of a gaussian convoluted with image. can take any value. 0= no smoothing. 

%===================================
% SVD analysis
%===================================
load('REF_Absorption_Spectra.mat')
clear REF_Spec LA 

Absorption=-log10(FImage);

% PARAMETERS:
% focus data on interesting part ofg spectrum, e.g. (500-700nm):
LambdaMin=475;      % min/max wavelengths of interest for the analysis
LambdaMax=700;

Ns=4;                % number of input spectra, including Hematoxylin, Eosin and two for background
dn=1;                % sub-sampling in X,Y  to reduce amount of computation by dn^2, if needed
Nr=50000;            % progress report rate
%...........................................................................

Ny=size(Absorption,1);
Nx=size(Absorption,2);
NL=size(Absorption,3);


LIdx=(lambda>LambdaMin).*(lambda<LambdaMax);
Lmin=find(LIdx, 1 ); 
Lmax=find(LIdx, 1, 'last' );             
L=lambda(Lmin:Lmax)';
NL=Lmax-Lmin+1;          


REF_Spec(:,1)=L;  %wavelength
% Actuall absorption curves of Eosin:
REF_Spec(:,2)=interp1(Eosin_Abs_Rel(:,1),Eosin_Abs_Rel(:,2),REF_Spec(:,1));
% Gaussian approximaltion to Eosin:
% WLm2=520;  Sig2=10; 
% REF_Spec(:,2)=exp(-(REF_Spec(:,1)-WLm2).^2/(2*Sig2^2));

% Actuall absorption curves of Hematoxylin:
REF_Spec(:,3)=interp1(Hematoxylin_Abs_Rel(:,1),Hematoxylin_Abs_Rel(:,2),REF_Spec(:,1));
% Gaussian approximaltion to Hematoxylin:
% WLm3=590;  Sig3=40; 
% REF_Spec(:,3)=exp(-(REF_Spec(:,1)-WLm3).^2/(2*Sig3^2));

% Additional curves to account for background and unknown materials
WLm4=700;  Sig4=50; 
REF_Spec(:,4)=exp(-(REF_Spec(:,1)-WLm4).^2/(2*Sig4^2));
WLm5=460;  Sig5=30; 
REF_Spec(:,5)=exp(-(REF_Spec(:,1)-WLm5).^2/(2*Sig5^2));

% figure(5); plot(REF_Spec(:,1),REF_Spec(:,2),REF_Spec(:,1),REF_Spec(:,3))
% .............................................

A=REF_Spec(:,2:Ns+1);               % This is the transfer matrix from fractions to combined spectrum
[U,S,V] = svd(A);                        % SVD decomposition
S1=[inv(S(1:Ns,1:Ns)) zeros(Ns,NL-Ns)];                % this is the "inverse" of S

UnMix=inv(V')*S1*inv(U);            % this is generalized inv(A), allowing to go from spectra to fractions

%========================================================================================
% PREFILTERING:

LA=zeros(Ny,Nx,NL);

% Optional smoothing of the data (each WL seperately) in (X,Y) before analysis :
for j=1:NL
    if Ssmooth>0
        LA(:,:,j)=SpatialFilter(Absorption(:,:,Lmin+j-1),Ssmooth);    
    else
        LA(:,:,j)=Absorption(:,:,Lmin+j-1);
    end
end

% Optional subsampling:
if dn~=1   
    LA=LA(1:dn:end,1:dn:end,:);
end

%========================================================================================
% Running SVD analysis per spectrum:

Fraction=zeros(Ny,Nx,Ns);
Err=zeros(Ny,Nx);

n=0;
Tstart=now;

for j=1:Ny
    for k=1:Nx
        n=n+1;
        S=squeeze(LA(j,k,:));
        fr=UnMix*S;
        Fraction(j,k,:)=fr;                    % Fraction of each component is the output of analysis
        Err(j,k)=std(S-A*fr)^2;
        
        if fix(n/Nr)==n/Nr                     % progress report
            pFinish=100*n/(Nx*Ny);
            tPass= (now-Tstart)*24*60;
            tLeft = tPass*(100-pFinish)/pFinish;
            if ProgRep=='t'
                % progress report by text:
                disp([num2str(pFinish) ' % ||  Time from start: ' num2str(tPass ) ' min. ||  Remaining Time: ' num2str(tLeft)  ' min.  '])
            elseif ProgRep=='g'
                % progress report by graph:
                figure(22); clf;
                C1=A(:,1)*fr(1);
                C2=A(:,2)*fr(2);
                C3=A(:,3)*fr(3);
                C4=A(:,4)*fr(4);
                plot(L,S,'or', L,A*fr,'-b',L,C1,'--c',L,C2,'--c',L,C3,'--c',L,C4,'--c')
                drawnow
                title(['Error (SSE): ' num2str( Err(Ny,Nx)) ' ||  Finished: ' num2str(pFinish) ' % ||  Time from start: ' num2str(tPass ) ' min. ||  Remaining Time: ' num2str(tLeft)  ' min.  ']);
                pause(0.1)
            end
        end
    end
end
end
function [IndAll,Spectra]=ShowSpec2(RGBHE,spec,lambda,n,Fraction)
% RGBHE - RGB Image for the original H&E transmission
% spec - data cube
% labmbda ...
% n - number of nearest neibors for averaging when displaying spectrum.

FRGB=FalseFractionRGB(Fraction,[0,0,3,0; 0,0.7,0,0; 1,0,0,1]);      % false RGB representing SVD analysis results (Fraction)

f=figure(6);
if isempty(f.CurrentAxes)               % if figure exist no need to redraw so magnif remains as before
    subplot(2,2,1); 
    imshow(RGBHE)  
    a=f.CurrentAxes;
    subplot(2,2,2);
    imshow(FRGB)
    b=f.CurrentAxes;
end

a.Position=[0.0    0.48    0.6    0.5];
b.Position=[0.5    0.48    0.6    0.5];

IndAll=[];
Spectra=[];

stop=0;

while stop~=1
    [x,y,button]=ginput(1);
    X=round(x);
    Y=round(y);
    Ind=[X Y];
    if button=='s' || button=='S'       % 's' for stopping
        stop=1; 
    elseif button=='p' || button=='P'   % 'p' may be used for pausing and changing magnification, hit Enter when done 
        pause
        % adjust magnif similar for both views:
        c=f.CurrentAxes;    % read imag limits from last touched subplot
        XL=c.XLim;          
        YL=c.YLim;
        a.XLim=XL;          % apply same limits to both plots
        a.YLim=YL;
        b.XLim=XL;
        b.YLim=YL;
    elseif button=='n' || button=='N'   % 'n' for changing the averaging area
        n=input('Averaging area size? : ');
    else                                % mouse-click - choosing a point a drawing spectrum
        subplot(2,2,3); 
        a3=f.CurrentAxes;
        a3.Position=[0.1300    0.0800    0.3347    0.3];
        if n==0
            Spectrum=squeeze(spec(Y,X,:));
            F1=Fraction(Y,X,1);
            FE=Fraction(Y,X,2);
            FH=Fraction(Y,X,3);
            F4=Fraction(Y,X,4);
        else
            Spectrum=squeeze(mean(spec(Y-n:Y+n,X-n:X+n,:),[1 2]));
            F1=mean(Fraction(Y-n:Y+n,X-n:X+n,1),[1 2]);
            FE=mean(Fraction(Y-n:Y+n,X-n:X+n,2),[1 2]);
            FH=mean(Fraction(Y-n:Y+n,X-n:X+n,3),[1 2]);           
            F4=mean(Fraction(Y-n:Y+n,X-n:X+n,4),[1 2]);           
        end
        plot(lambda,Spectrum)
        axis([400 800 0 1.2])
        title(['X:  ' num2str(X) ' , Y:  ' num2str(Y)])
        IndAll=[IndAll; Ind];
        Spectra=[Spectra; Spectrum];
        
        subplot(2,2,4)
        a4=f.CurrentAxes;
        a4.Position=[0.5703    0.1100    0.3347    0.3412];

        bar(1:4,[F1 FE FH F4])
        axis([0 5 -.2 1.2])
        title(['F1:  ' num2str(F1,'%.3f')   '      , Eosin:  ' num2str(FE,'%.3f') '       ,  Hematoxylin:  ' num2str(FH,'%.3f') '      , F4:  ' num2str(F4,'%.3f') ])
        text(1,1,['n= ' num2str(n)])
  
    end
   
end

end
function Af=SpatialFilter(A,Nf)

% Nf - filter half-size:
x=-Nf:Nf;
y=x';
X=ones(2*Nf+1,1)*x;
Y=y*ones(1,2*Nf+1);
F=exp(-2*(X.^2+Y.^2)/Nf^2);        % Gaussian filter definition
F=F/sum(sum(F));         

Af=conv2(A,F,'same');        % applying filter to Mask
end
function FRGB=FalseFractionRGB(Fraction,W)

% W = weights for RGB false coloring
% e.g. W=[0,0,2,0; 0,0.5,0,0; 1,0,0,1]


FRGB(:,:,1)=W(1,1)*Fraction(:,:,1)+W(1,2)*Fraction(:,:,2)+W(1,3)*Fraction(:,:,3)+W(1,4)*Fraction(:,:,4);      % R
FRGB(:,:,2)=W(2,1)*Fraction(:,:,1)+W(2,2)*Fraction(:,:,2)+W(2,3)*Fraction(:,:,3)+W(2,4)*Fraction(:,:,4);      % G
FRGB(:,:,3)=W(3,1)*Fraction(:,:,1)+W(3,2)*Fraction(:,:,2)+W(3,3)*Fraction(:,:,3)+W(3,4)*Fraction(:,:,4);      % B
end
function XRGB=Spec2RGB(spec,lambda)

R=find(lambda>600);
G=find((lambda>500).* (lambda<600));
B=find(lambda<500);
Dim=size(spec);

XRGB(:,:,1)=sum(spec(:,:,R),3);
Map=squeeze(XRGB(:,:,1));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,1)=XRGB(:,:,1)./prctile(Sig,98);

XRGB(:,:,2)=sum(spec(:,:,G),3);
Map=squeeze(XRGB(:,:,2));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,2)=XRGB(:,:,2)./prctile(Sig,98);

XRGB(:,:,3)=sum(spec(:,:,B),3);
Map=squeeze(XRGB(:,:,3));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,3)=XRGB(:,:,3)./prctile(Sig,98);
end
%}
%Treshold + histogram applied for h and for e
%{
h = imread('current_work_h.png');
e = imread('current_work_e.png');
q = e<150;
imtool(q);
histogram(h(:));
%}
%-----------------------------------------------------------------

%SMOOTHING WL x40

%{
%some general filters
weighted_filter = [1 2 1; 2 4 2; 1 2 1]/16;
horiz_filter = [0 0 0; 1 1 1; 0 0 0];
vert_filter = [0 1 0; 0 1 0; 0 1 0];
ones_filt = [ 1 1 1; 1 1 1; 1 1 1];
mean_filter_5 = ones(5,5)/25;
%meadian_filter = medfilt2(img);
gaussian = fspecial('gaussian');
laplacian = fspecial('laplacian');
log = 
%remind step3 remove noise
%{
bw = bwareaopen(bw,30);
%}

g = load('x40_r1_r2_spec.mat');
lambda=g.lambda; 
lambda = lambda(1:28);
wl_mat =g.spec;
chosen_wl=1;
grayscale_out = read_unique_wl(wl_mat , lambda, chosen_wl);
%imtool(grayscale_out);

%smoothing
%{
smoothed = imfilter(grayscale_out,weighted_filter, 'same');
%or conv2 or filter2 for conv
imtool(uint8(255 * mat2gray(smoothed)));
%}


%  same smoothing to all WL , and  SAVE the new spectrum + its lambda
%{
[a,b] = size( grayscale_out);
smoothed = imfilter(grayscale_out,weighted_filter, 'same');
new_mat = smoothed;

for i =2:28
chosen_wl=i;
grayscale_out = read_unique_wl(wl_mat , lambda, chosen_wl);
smoothed = imfilter(grayscale_out,weighted_filter, 'same');
new_mat = cat(3,new_mat, smoothed);
end
size(new_mat)
smoothed_x40  = cast(new_mat, 'double');
rgb = spectoRGB(smoothed_x40, lambda, 0,'www');
white_rgb = whitening_rgb (rgb,412,1574,567,1723,3760, 155,3946,281,0,'www');

name = "x40_r1_r2_smoothed_weighted_filter";
save(name,'smoothed_x40','lambda','rgb', 'white_rgb');
%}

%}

%-----------------------------------------------------------------

% Watershed technique
%{
%https://fr.mathworks.com/help/images/marker-controlled-watershed-segmentation.html


rgb = imread('C161VA_x20_X700Y600_Whitened_2.png');
rgb = imcrop(rgb,[1500 600  (2000-1500) (950-600)]);
I = rgb2gray(rgb);
I2 = imcomplement(I);
%imtool(I);
%case: kmeans instead of rgb
%{
kmeans_20 = imread('C161VA_x20_X700Y600_20_kmeans.png');
kmeans_20 = imcrop(kmeans_20,[1500 600  (2000-1500) (950-600)]);
I = kmeans_20;
I2= kmeans_20;
%}



%gradient magnitude as segmentation function
gmag = imgradient(I);
imshow(gmag,[])
title('Gradient Magnitude')
%case: using edge instead of gradients.
%{
hne= imread('c161_HnE.png');
hne = imcrop(hne,[1500 600  (2000-1500) (950-600)]);

hne=rgb2gray(hne);
hnedge = edge(hne);

c_40 = imread('C161VA_x20_X700Y600_20_kmeans_c_40_0_3000_on_empty.png');
c_40 = imcrop(c_40,[1500 600  (2000-1500) (950-600)]);
c40edge=edge(c_40);

c_242= imread('C161VA_x20_X700Y600_20_kmeans_c_242_20_3000_on_k13_20_1000.png');
c_242 = imcrop(c_242,[1500 600  (2000-1500) (950-600)]);
c242edge=edge(c_242);

%masked_img=bw_mask_on_bw(hnedge, edge(c_40));
masked_img=bw_mask_on_bw(hnedge, c40edge);
masked_img=bw_mask_on_bw(masked_img, c242edge);
imtool(masked_img);

gmag = uint8(masked_img);
%}


%mark foreground objct with morphological techniques:  "opening-by-reconstruction" and "closing-by-reconstruction" to "clean" up the image. 
%basic opening
%{
%Opening = erosion then dilation
se = strel('disk',4);
Io = imopen(I2,se);
imshow(Io)
title('Opening')
%}

% opening-by-reconstruction 
se = strel('disk',4);
Ie = imerode(I2,se);
Iobr = imreconstruct(Ie,I2);
imshow(Iobr);
title('Opening-by-Reconstruction')


%basic closing
%{
%Following the opening with a closing can remove the dark spots and stem marks. Compare a regular morphological closing with a closing-by-reconstruction. First try imclose:
Ioc = imclose(Io,se);
imshow(Ioc)
title('Opening-Closing')
%}


% closing-by-reconstruction 
%Now use imdilate followed by imreconstruct. Notice you must complement the image inputs and output of imreconstruct.

Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
imshow(Iobrcbr)
title('Opening-Closing by Reconstruction')
%As you can see by comparing Iobrcbr with Ioc, reconstruction-based opening and closing are more effective than standard opening and closing at removing small blemishes without affecting the overall shapes of the objects.


%Calculate the regional maxima of Iobrcbr to obtain good foreground markers.
fgm = imregionalmax(Iobrcbr);
imshow(fgm)
title('Regional Maxima of Opening-Closing by Reconstruction')


%To help interpret the result, superimpose the foreground marker image on the original image.
I3 = labeloverlay(I2,fgm);
imshow(I3)
title('Regional Maxima Superimposed on Original Image')


%not wanted process
%{
%Notice that some of the mostly-occluded and shadowed objects are not marked,
%which means that these objects will not be segmented properly in the end result. 
%Also, the foreground markers in some objects go right up to the objects' edge.
%That means you should clean the edges of the marker blobs and then shrink them a bit.
%You can do this by a closing followed by an erosion.

se2 = strel(ones(5,5));
fgm2 = imclose(fgm,se2);
fgm3 = imerode(fgm2,se2);

%This procedure tends to leave some stray isolated pixels that must be removed. 
%You can do this using bwareaopen, which removes all blobs that have fewer than a certain number of pixels.

fgm4 = bwareaopen(fgm3,5);
I4 = labeloverlay(I2,fgm4);
imshow(I4)
title('Modified Regional Maxima Superimposed on Original Image')

%}


%Step 4: Compute Background Markers
%Now we mark the background. In the cleaned-up image, Iobrcbr, the dark pixels belong to the background, so start with a thresholding operation.

bw = imbinarize(Iobrcbr);
imshow(bw)
title('Thresholded Opening-Closing by Reconstruction')

%The background pixels are in black, but ideally we don't want the background markers to be too close to the edges of the objects we are trying to segment. We'll "thin" the background by computing the "skeleton by influence zones", or SKIZ, of the foreground of bw. This can be done by computing the watershed transform of the distance transform of bw, and then looking for the watershed ridge lines (DL == 0) of the result.
%{
D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;
imshow(bgm)
title('Watershed Ridge Lines')
%}
%case c40_242 used for skeleton

c40_242 = imread('C161VA_x20_X700Y600_20_kmeans_c_40_242_15_15000_on_empty.png');
c40_242 = imcrop(c40_242,[1500 600  (2000-1500) (950-600)]);
%imtool(c40_242);

D = bwdist(c40_242);
DL = watershed(D);
bgm = DL == 0;
imshow(bgm)
title('Watershed Ridge Lines')



%Step 5: Compute the Watershed Transform of the Segmentation Function.
%The function imimposemin can be used to modify an image so that it has regional minima only in certain desired locations. Here you can use imimposemin to modify the gradient magnitude image so that its only regional minima occur at foreground and background marker pixels.

gmag2 = imimposemin(gmag, bgm | fgm);

%Finally, compute the watershed-based segmentation.
L = watershed(gmag2);

%Step 6: Visualize the Result
%One visualization technique is to superimpose the foreground markers, background markers, and segmented object boundaries on the original image. You can use dilation as needed to make certain aspects, such as the object boundaries, more visible. Object boundaries are located where L == 0. The binary foreground and background markers are scaled to different integer values so that they are assigned different labels.
labels = imdilate(L==0,ones(3,3)) + 2*bgm + 3*fgm;
I4 = labeloverlay(I,labels);
imshow(I4)
title('Markers and Object Boundaries Superimposed on Original Image')

%This visualization illustrates how the locations of the foreground and background markers affect the result. In a couple of locations, partially occluded darker objects were merged with their brighter neighbor objects because the occluded objects did not have foreground markers.

%Another useful visualization technique is to display the label matrix as a color image. Label matrices, such as those produced by watershed and bwlabel, can be converted to truecolor images for visualization purposes by using label2rgb.
Lrgb = label2rgb(L,'jet','w','shuffle');
imtool(Lrgb)
title('Colored Watershed Label Matrix')

%You can use transparency to superimpose this pseudo-color label matrix on top of the original intensity image.
figure
imshow(I)
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title('Colored Labels Superimposed Transparently on Original Image')

%}
%-----------------------------------------------------------------
%R1 R2 merge x
%{
R1 and R2 loading
%{
mymat1 = load('C161VA_X40_X600Y0_R1.mat');
rgb1 = mymat1.RGB;
imtool(rgb1);

mymat2 = load('C161VA_X40_X600Y0_R2.mat');
rgb2 = mymat2.RGB;
imtool(rgb2);
%}
%{
%{
%mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/3.rawstack.mat');

mymat = load('C161VA_X40_X600Y0_R2.mat');

input_spectral_img = mymat.spec(:,:,:);
number_of_clusters = 20;
%1 if want to save img
save_kmean_img=0;
name_of_input_img = 'x40';
%}

%{
noiseless = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/c161va_without_noise.mat');
input_spectral_img = noiseless.FImage;
number_of_clusters = 20;
save_kmean_img=0;
name_of_input_img = 'noiseless';
%}

kmean_mat_output=kmeans_algo(input_spectral_img, number_of_clusters, save_kmean_img, name_of_input_img);
imtool(kmean_mat_output);

%}
%}
%-----------------------------------------------------------------

%merging edges of hne , 40 and 242, half thicken them , then apply filter
%detection on them in order to eventually get cells segm.
%{
%{
%{  
                       %!
hne= imread('c161_HnE.png');
hne = imcrop(hne,[1500 600  (2000-1500) (950-600)]);

hne=rgb2gray(hne);
hnedge = edge(hne);

c_40 = imread('C161VA_x20_X700Y600_20_kmeans_c_40_0_3000_on_empty.png');
c_40 = imcrop(c_40,[1500 600  (2000-1500) (950-600)]);
c40edge=edge(c_40);

c_242= imread('C161VA_x20_X700Y600_20_kmeans_c_242_20_3000_on_k13_20_1000.png');
c_242 = imcrop(c_242,[1500 600  (2000-1500) (950-600)]);
c242edge=edge(c_242);

%masked_img=bw_mask_on_bw(hnedge, edge(c_40));
masked_img=bw_mask_on_bw(hnedge, c40edge);
masked_img=bw_mask_on_bw(masked_img, c242edge);
%imtool(masked_img);


b = uint8(~masked_img);
imtool(mat2gray(b));
%}

%edge thicker
%{
ba = uint8(masked_img);

b2 = bwmorph(b,'thin',2); %thicken
imtool(mat2gray(b2));

%b3 = imdilate(ba, strel('disk',1));
%imtool(mat2gray(b3));
%}


%half thiker            !

[k,l]=size(b);
c =zeros(k,l);
for i=2:k-1
    for j=2:l-1
        if ~b(i,j)==1
           % c(i-1, j-1)=1;
      %      c(i-1, j)=1;
           % c(i-1, j+1)=1;
            c(i, j-1)=1;
            c(i, j)=1;
       %     c(i, j+1)=1;  
           % c(i+1, j-1)=1;
            c(i+1, j)=1;
          %  c(i+1, j+1)=1;
        end
    end
end
        
imtool(mat2gray(c));
%}
% c_smooth = conv(c, ones(1,5)/5, 'same');
%bluring 
%{
windowSize = 2; % Bigger for more blurring.
kernel = ones(windowSize) / windowSize^2;
output = imfilter(c, kernel);
imtool(output);
%}
%if_i                  !
 %{
 % SELECT working_img : which kmeans to work on . 
%{
%temp = imread('10_kmeans_of_mat_P055A_X2000Y0_R1.png');
%crop_k =mat2gray(temp);

%first img loads
temp = imread('10K_means_out.png');
cropk= imcrop(temp,[520 400  (819-520) (699-400)]);
%imtool(cropk);
%}

%reference cell
wkn_img = ~c;
%wkn_img=imcrop(wkn_img,[1500 600  (2000-1500) (950-600)]);
name = 'c161va_edges_40_242_hne';

% select  clusters to filter 
cell_area = [1];
cell_area_name = 'edges';
% select min - max allowed size
Max_allowed_cc_size = 4000;
Min_allowed_cc_size = 15;
%select background
%{
%kernel_out_8_300 = imread('kernel_out_8_300.png');
%support_name= 'kernel_8_300';
%support=mat2gray(kernel_out_8_300);

%m_k_8_300 = imread('kernel_8_300_on_empty_filtered_img.png');
%support_name= 'm_k_8_300';
%support=mat2gray(m_k_8_300);

%membrane_kernel_8_100 = imread('kernel_out_8_100.png');
%support =mat2gray(membrane_kernel_8_100);
%support_name= 'membrane_kernel_8_100';


[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;


background = mat2gray(imread('C161VA_x20_X700Y600_20_kmeans_k_13_20_1000_on_empty.png'));
bg_name= 'k13_20_1000';

%}
[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;
%last param
bg_reduction = 0;
save_filtered_img = 0;
final_name = append(name, '_', cell_area_name,'_',int2str(Min_allowed_cc_size),'_',int2str(Max_allowed_cc_size),'_on_',bg_name);

f_i = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction,4);

 mask_kernel= f_i;
 imtool(mask_kernel)
 %}
 
%}


%-----------------------------------------------------------------

%from kmeans clusters to circles on rgb , using app circle finder
%imtool(background_full);
%imtool(background);
%{
tissue_name ='x40_r1_r2_spec';
struct_path = '/Users/lior/Desktop/Image & Analysis /' +  string(tissue_name) + '.mat';
data = load(struct_path);


crop_area_1 = [500 1000  (1000-500) (1500-1000)];
crop_area_2 = [1750 580  (2150-1750) (820-580)];
crop_area_3 = [ 1500 1000  (2000-1500) (1500-1000)];
crop_area_4 = [ 3400 1000  (3800-3400) (1400-1000)];
img4 = imcrop(data.rgb, crop_area_4);


%qwert = insertShape(background_full,'filledrectangle', crop_area_1, 'Color', 'blue','Opacity',.3);
%qwert = insertShape(qwert,'filledrectangle', crop_area_2, 'Color', 'blue','Opacity',.3);
%qwert = insertShape(qwert,'filledrectangle', crop_area_3, 'Color', 'blue','Opacity',.3);
%qwert = insertShape(qwert,'filledrectangle', crop_area_4, 'Color', 'blue','Opacity',.3);

%}


%select pxl manually : prototype of fct manual_cell_labelization
%{
tissue_name = 'x40_r1_r2';
focused_area_name = 'area_4';
x=0;
b=false;
cell_in_area_4 = [];
while x<100 && b==false
    imshow(img4);
    roi = drawellipse(); % select manually gt
    x = x+1;
    save('/Users/lior/Desktop/Image & Analysis /cells_label/'+ string(tissue_name) +'/'+ string(focused_area_name) + '/cell_' + string(num2str(x)) + '.mat', 'bw');

    %cell_in_area_4 = [cell_in_area_4 ; roi.Vertices];
    %cell_in_area_4 = [cell_in_area_4 ; -1,-1];
    

end
%save('cell_in_area_4_'+string(tissue_name)+'.mat' ,'cell_in_area_4');


bw = createMask(roi);
q=1;
%}



%suite process after selection and saving labbelled cell . 
%{

%select first roi by stopping via (-1,-1) rows.
labeled_cell_1 = [];
[s,~]  =  size(cell_in_area_4);
i=1;
coord = round(cell_in_area_4(i,:));
while coord(1) ~= -1 && coord(2) ~= -1 && i<= s
    labeled_cell_1 = [labeled_cell_1; coord(1),coord(2)];
    i=i+1;
    coord  = round(cell_in_area_4(i,:));
end
[size_lab_cell_1,~] = size(labeled_cell_1);




labelled_area_4 = img4;
for j=1:size_lab_cell_1
    labelled_area_4(labeled_cell_1(j,1),labeled_cell_1(j,2),:) = [0.100,0.100,0.100];
end
imshow(labelled_area_4);
%}

%other & previous tries for selectied manually circles in img
%{
f1_mserR_cluster1_a1.PixelList{1, 1}  % pred location of pixels.
imtool(distance_map_list_a1);
imtool(distance_map_list_a2);
imtool(distance_map_list_a3);
%img = imread('/Users/lior/Documents/MATLAB/ex4_image_processing/balls1.tiff');
%while true
img = distance_map_list_a1;
imshow(img);
roi = drawcircle('Color','r');

% img = add_circles_to_img(img, roi.Center, roi.Radius);
imtool(img);
%end
x=0;




% Create a logical image of an ellipse with specified
% semi-major and semi-minor axes, center, and image size.
% First create the image.
imageSizeX = 501;
imageSizeY = 501;
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the ellipse in the image.
centerX = 253;
centerY = 303;
radiusX = 250;
radiusY = 150;
ellipsePixels = (rowsInImage - centerY).^2 ./ radiusY^2 ...
    + (columnsInImage - centerX).^2 ./ radiusX^2 <= 1;
% ellipsePixels is a 2D "logical" array.
% Now, display it.
image(ellipsePixels) ;
colormap([0 0 0; 1 1 1]);
title('Binary image of a ellipse', 'FontSize', 20);


%{
centers_list = [
320,250, 61
341,129, 50
418, 156 ,43
308 , 19 , 18


];
radius_list = [61, 50];
%}
centx=100;
centy=100;
r =40;
%[xgrid, ygrid] = meshgrid(1:size(img,2), 1:size(img,1));
%mask = ((xgrid-x).^2 + (ygrid-y).^2) <= r.^2;
%values = img(mask);

%{
x_list = [];
y_list = [];
 if true
    theta = 0 : (2 * pi / 10000) : (2 * pi);
    for elt= 1:length(theta):
        if floor(cos(elt)) == cos(elt) && floor(sin(elt)) == sin(elt):
             x_list = [x_list ; r * cos(theta) + centx];
        end
        if floor(sin(elt)) == sin(elt):
             y_list = [y_list ; r * sin(theta) + centy];
        end
    end
    pline_x =uint8( r * cos(theta) + centx);
    pline_y = round(r * sin(theta) + centy,0);
end



%}



%imtool(img);
%h = images.roi.Circle(gca,'Center',[1000 1000],'Radius',500);

%}
% end of tries for cell gt labelization
%--------------------------------------------------------

%{

%edge detect from f_i clusters of kmeans.
edge_factor = .013; %since edge detect is onkmeans, edge_factor param is useless
c_40_242 = imread('C161VA_x20_X700Y600_20_kmeans_c_40_242_15_150000000_on_empty.png');
c_crop = imcrop(c_40_242,[1500 600  (2000-1500) (950-600)]);

%c_crop = imcrop(c_40_242,[1700 850  (1815-1700) (930-850)]);

mask_kernel = edge(c_crop,'sobel',edge_factor);
imtool(mask_kernel);
%convex object in img
%{
k  = bwconvhull(mask_kernel,'objects');
imtool(k);
%}
 %if needed f_i
 %{
 % SELECT working_img : which kmeans to work on . 
%{
%temp = imread('10_kmeans_of_mat_P055A_X2000Y0_R1.png');
%crop_k =mat2gray(temp);

%first img loads
temp = imread('10K_means_out.png');
cropk= imcrop(temp,[520 400  (819-520) (699-400)]);
%imtool(cropk);
%}

%reference cell
wkn_img = edges;
%wkn_img=imcrop(wkn_img,[1500 600  (2000-1500) (950-600)]);
name = 'c161va_edge_sobel_on_cancerousref1';

% select  clusters to filter 
cell_area = [1];
cell_area_name = 'edges';
% select min - max allowed size
Max_allowed_cc_size = 1000000000;
Min_allowed_cc_size = 0;
%select background
%{
%kernel_out_8_300 = imread('kernel_out_8_300.png');
%support_name= 'kernel_8_300';
%support=mat2gray(kernel_out_8_300);

%m_k_8_300 = imread('kernel_8_300_on_empty_filtered_img.png');
%support_name= 'm_k_8_300';
%support=mat2gray(m_k_8_300);

%membrane_kernel_8_100 = imread('kernel_out_8_100.png');
%support =mat2gray(membrane_kernel_8_100);
%support_name= 'membrane_kernel_8_100';


[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;


background = mat2gray(imread('C161VA_x20_X700Y600_20_kmeans_k_13_20_1000_on_empty.png'));
bg_name= 'k13_20_1000';

%}
[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;
%last param
bg_reduction = 0;
save_filtered_img = 0;
final_name = append(name, '_', cell_area_name,'_',int2str(Min_allowed_cc_size),'_',int2str(Max_allowed_cc_size),'_on_',bg_name);

f_i = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction);

 mask_kernel= f_i;
 %}

% edge mask over rgb img
whit = imread('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/C161VA_x20_X700Y600_Whitened_2.png');
viz = imcrop(whit,[1500 600  (2000-1500) (950-600)]);
%viz =imcrop(imread('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/c161_HnE.png'),[1500 600  (2000-1500) (950-600)]);

masked_img=red_mask_on_img(viz, mask_kernel);
imtool(masked_img);

%put circles founded by circle finder on img
img = imread('C161VA_x20_X700Y600_Whitened_2.png');
c_crop = imcrop(img,[1500 600  (2000-1500) (950-600)]);
centers_list = tcs; % tcs and rcs a got from circlefinder app . params to bet set . 
radii_list = trs;

img_with_circles= add_circles_to_img(c_crop, centers_list, radii_list);
imtool(img_with_circles);
%}

%-----------------------------------------------------------------
%temporary strip aliner:

%{
r1 = load('C161VA_X40_X600Y0_R1.mat');
r1spec = r1.spec;
r2 =load('C161VA_X40_X600Y0_R2.mat');
r2spec= r2.spec;
[q,w,e] = size(r1spec);
[a,z,x] = size(r2spec);
%166 : manually chosen strip aligner: 166 last lines.
r1s =r1spec(1:q-83,1:w-34,:);
r2s = r2spec(84:a, 35:z , :);
lambda=r2.lambda; 

imtool(Spec2RGB(r1spec , lambda));
imtool(Spec2RGB(r2spec, lambda));
spec = cat(1,r1s, r2s);
%size(x40_merged)
rgb = Spec2RGB(spec,lambda);
imtool(rgb);
save( 'x40_r1_r2_spec','spec' ,'lambda','rgb', '-v7.3');

%}
%WRONG  strip aliner
%{
%{

%not good aligned
r1 = load('C161VA_X40_X600Y0_R1.mat');
r1spec = r1.spec;
r2 =load('C161VA_X40_X600Y0_R2.mat');
r2spec= r2.spec;
ch=r2.lambda; 
%imtool(Spec2RGB(r2spec , ch));
%imtool(Spec2RGB(r3spec, ch));

%34
CombImage=ImageAligner(r1spec,r2spec,34);
imtool(Spec2RGB(CombImage,ch));
%}

function CombImage=ImageAligner(Top,Bot,Ch)
% Stitch together two image stripes

global lambda

% define size of target for alignment in bottom of Top
TargXMinF=0.3;      % part of Top width
TargXMaxF=0.7;      % part of Top width
TargY=50;          % heigh, in pixels

% Top,Bot - input matrices (2 or 3 dimensions) according to their intended positions
% Ch -  the number of the channel (3rd dim) according to which stiching is done. 
    %   If the image is 2 dimensional (no color or wavelenght) use Ch=1
% 
XL=size(Top,2);
YL=size(Top,1);

if ndims(Top)>2
    TopI=Top(:,:,Ch);
    BotI=Bot(:,:,Ch);
    Nc=size(Top,3);     % number of channels (wavelengths/colors)
else
    TopI=Top;
    BotI=Bot;
    Nc=1;
end

figure(10); clf;

subplot(2,1,1); 
imshow(TopI/max(max(TopI)));
title('Top')
hold on
subplot(2,1,2); 
imshow(BotI/max(max(BotI)));
title('Bot')
hold on
% sgtitle('Mark top-left & bot-right of a target in TOP image')

% [xc,yc]=ginput(2);
xc=[round(TargXMinF*XL) round(TargXMaxF*XL)];
yc=[YL-TargY YL];

Target=TopI(floor(yc(1)):floor(yc(2)),floor(xc(1)):floor(xc(2)));
Target=Target-mean(mean(Target));       % removing offset cancels the possibility to get false correlation from large bright areas

nBot=BotI-mean(mean(BotI));
nTop=TopI-mean(mean(TopI));

% corr1=xcorr2(nTop,Target);
% [Y1,X1]=find(corr1==max(max(corr1)))
X1=round(xc(2));
Y1=round(yc(2));
corr2=xcorr2(nBot(1:500,:),Target);
[Y2,X2]=find(corr2==max(max(corr2)));

% figure(2);
subplot(2,1,1);
% imshow(TopI/max(max(TopI)));
% imshow(corr1/max(max(corr1)));
% hold on
plot(X1,Y1,'r+')
% title('Top')
subplot(2,1,2);
% imshow(BotI/max(max(BotI)));
% imshow(corr2/max(max(corr2)));
% hold on
plot(X2,Y2,'r+')
% title('Bot')

dX=X1-X2;
dY=Y1-Y2;

W=size(Top,1);
Lt=size(Top,2);
Lb=size(Bot,2);

if ndims(Top)>2
    TmpTop=ones(dY,Lb,Nc);
    if dX>=0
        TmpTop(: , (1:(Lb-dX)),:)= Top( (1:dY),(dX+1:Lb),: );
    else
        TmpTop(: , (-dX+1:Lb),:)= Top( (1:dY),(1:Lb+dX),: );
    end
    CombImage = permute([ permute(TmpTop,[2 1 3]) permute(Bot,[2 1 3]) ],[2 1 3]);
else
    if dX>=0
        TmpTop(: , (1:(Lb-dX)))= Top( (1:dY),(dX+1:Lb) );
    else
        TmpTop(: , (-dX+1:Lb))= Top( (1:dY),(1:Lb+dX) );
    end
    CombImage = permute([ permute(TmpTop,[2 1]) permute(Bot,[2 1]) ],[2 1]);
end

figure(11); clf;

if ndims(Top)==2                 % B/W image
    imshow(CombImage)
elseif size(Top,3)==3           % RGB image
    imshow(CombImage)
else                            % full spectrum-image
    imshow(Spec2RGB(CombImage,lambda))
end


end

%}

%{
healthy = imread('C161VA_x20_X700Y600_20_kmeans_k_13_20_1000_on_empty.png');
cancerous = imread('C161VA_x20_X700Y600_cancerous_zone.png');
both = healthy + cancerous;
imtool(both);

%}

%-----------------------------------------------------------------

% test on noise removing boaz
%{
%{
mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/EmptyRef.mat');
Empty = mymat.spec;
lambda = mymat.lambda;
EmptyRef = FixedPatternLearn(Empty,14,3);
size(EmptyRef)
%imtool(EmptyRef);
%{
img = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/C161VA_x20_X700Y600.mat');
immg = img.spec;
RGBHE =  imread('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/C161VA_x20_X700Y600_SpecToRgb.png');
FImage = FixedPatternCorrect(immg, lambda,EmptyRef, 1, 0);
%imtool(FImage(:,:,4)/max(max(FImage(:,:,4))));
%imtool(FImage);
%save('c161va_without_noise.png', 'FImage');

[Fraction,Err]=SVDAnalysisHnE(FImage,lambda,'g',1);
%FRGB=FalseFractionRGB(Fraction,[0,0,3,0; 0,0.7,0,0; 1,0,0,1]);      % false RGB representing SVD analysis results (Fraction)
%imtool(FRGB);
%imwrite(FRGB, 'c161_HnE.png');
[IndAll,Spectra]=ShowSpec2(RGBHE,immg,lambda,0,Fraction);
%}
%}
function FImage=FixedPatternCorrect(Image,lambda,EmptyRef,SmoothLength,WB)
%{
% Image - for correction
% lambda - vector of wavelengths
% EmptyRef - a reference matrix produced from an Empty screen (or part of
%       image) using FixedPatternLearn.m
% SmoothLength: e.g. =1 smoothing along Y before applying calibration with EmptyRef
% WB: 1 - perform White balance
NwlMax=size(Image,3);
MeasSkipLength=size(EmptyRef,2);       % the skip length during the measurement 
FImage=zeros(size(Image));


% Optional smoothing along Y axis - improve SNR but may hurt resolution:
if SmoothLength>0
    SmoothVector=[(1:1:SmoothLength) (SmoothLength:-1:1)];  
    SmoothVector=SmoothVector'/sum(SmoothVector);
    SmthImage=zeros(size(Image));
    for Nwl=1:NwlMax
         SmthImage(:,:,Nwl)=conv2(Image(:,:,Nwl),SmoothVector,'same');
    end
else
    SmthImage=Image;
end

% normalizing the whole image accoring to reference provided (EmptyRef):   

for Nwl=1:NwlMax     % loop over wavelengths
           
    for j=1:size(Image,2)
        k=mod(j,MeasSkipLength);        % finding the correct normalization vector
        if k==0; k=MeasSkipLength; end
        FImage(:,j,Nwl)=SmthImage(:,j,Nwl)./EmptyRef(:,k,Nwl);
    end
    
end
%{
RGBHE=Spec2RGB(FImage,lambda);

if WB==1
    FRGB1=RGBHE;
    figure
    subplot(1,2,1);
    imshow(FRGB1)
    h=gca;
    h.Title.String='mark white area';
    rect=round(getrect);

    for c=1:3
        WhiteBal(c)=mean(mean(FRGB1(rect(2)+[1:rect(4)],rect(1)+[1:rect(3)],c)));
        RGBHE(:,:,c)=FRGB1(:,:,c)*(1/WhiteBal(c));
    end
    
    subplot(1,2,2);
    imshow(RGBHE)
end
%}

end
%}
function EmptyRef=FixedPatternLearn(RefImage,MeasSkipLength,YSmoothLength)
%{
% MeasSkipLength=14;       % skip length along X during the measurement 
% YSmoothLength=12;        % half length for smoothing along Y

NwlMax=size(RefImage,3);

SmoothVector=[(1:1:YSmoothLength) (YSmoothLength:-1:1)];  
SmoothVector=SmoothVector'/sum(SmoothVector);

Simage=zeros(size(RefImage));
EmptyRef=zeros(size(RefImage,1),MeasSkipLength,NwlMax);

% finding normalization matrix:

for Nwl=1:NwlMax     % loop over wavelengths
    
    % average along Y axis - improve SNR if no info in the figure (e.g. empty slide or empty part of figure)
    Simage(:,:,Nwl)=conv2(RefImage(:,:,Nwl),SmoothVector,'same');
    
    % average in jumps equivalent to X skipping.  notice this is done for full
    % vectors hence for each raw along the X scanning direction data is
    % collected separately in order to deal with nonuniformity along y axis.
    for StartCol=1:MeasSkipLength
        EmptyRef(:,StartCol,Nwl)=mean(Simage(:,(StartCol:MeasSkipLength:size(RefImage,2)),Nwl),2);
    end
       
end
end

%}
function [Fraction,Err]=SVDAnalysisHnE(FImage,lambda,ProgRep,Ssmooth)

% ProgRep='t';         % type of progress report:  t=text, g=graph, other= none
% Ssmooth=1;           % Smoothing on WL basis befor analysis, the values is the sigma of a gaussian convoluted with image. can take any value. 0= no smoothing. 

%===================================
% SVD analysis
%===================================
load('REF_Absorption_Spectra.mat')
clear REF_Spec LA 

Absorption=-log10(FImage);

% PARAMETERS:
% focus data on interesting part ofg spectrum, e.g. (500-700nm):
LambdaMin=475;      % min/max wavelengths of interest for the analysis
LambdaMax=700;

Ns=4;                % number of input spectra, including Hematoxylin, Eosin and two for background
dn=1;                % sub-sampling in X,Y  to reduce amount of computation by dn^2, if needed
Nr=50000;            % progress report rate
%...........................................................................

Ny=size(Absorption,1);
Nx=size(Absorption,2);
NL=size(Absorption,3);


LIdx=(lambda>LambdaMin).*(lambda<LambdaMax);
Lmin=find(LIdx, 1 ); 
Lmax=find(LIdx, 1, 'last' );             
L=lambda(Lmin:Lmax)';
NL=Lmax-Lmin+1;          


REF_Spec(:,1)=L;  %wavelength
% Actuall absorption curves of Eosin:
REF_Spec(:,2)=interp1(Eosin_Abs_Rel(:,1),Eosin_Abs_Rel(:,2),REF_Spec(:,1));
% Gaussian approximaltion to Eosin:
% WLm2=520;  Sig2=10; 
% REF_Spec(:,2)=exp(-(REF_Spec(:,1)-WLm2).^2/(2*Sig2^2));

% Actuall absorption curves of Hematoxylin:
REF_Spec(:,3)=interp1(Hematoxylin_Abs_Rel(:,1),Hematoxylin_Abs_Rel(:,2),REF_Spec(:,1));
% Gaussian approximaltion to Hematoxylin:
% WLm3=590;  Sig3=40; 
% REF_Spec(:,3)=exp(-(REF_Spec(:,1)-WLm3).^2/(2*Sig3^2));

% Additional curves to account for background and unknown materials
WLm4=700;  Sig4=50; 
REF_Spec(:,4)=exp(-(REF_Spec(:,1)-WLm4).^2/(2*Sig4^2));
WLm5=460;  Sig5=30; 
REF_Spec(:,5)=exp(-(REF_Spec(:,1)-WLm5).^2/(2*Sig5^2));

% figure(5); plot(REF_Spec(:,1),REF_Spec(:,2),REF_Spec(:,1),REF_Spec(:,3))
% .............................................

A=REF_Spec(:,2:Ns+1);               % This is the transfer matrix from fractions to combined spectrum
[U,S,V] = svd(A);                        % SVD decomposition
S1=[inv(S(1:Ns,1:Ns)) zeros(Ns,NL-Ns)];                % this is the "inverse" of S

UnMix=inv(V')*S1*inv(U);            % this is generalized inv(A), allowing to go from spectra to fractions

%========================================================================================
% PREFILTERING:

LA=zeros(Ny,Nx,NL);

% Optional smoothing of the data (each WL seperately) in (X,Y) before analysis :
for j=1:NL
    if Ssmooth>0
        LA(:,:,j)=SpatialFilter(Absorption(:,:,Lmin+j-1),Ssmooth);    
    else
        LA(:,:,j)=Absorption(:,:,Lmin+j-1);
    end
end

% Optional subsampling:
if dn~=1   
    LA=LA(1:dn:end,1:dn:end,:);
end

%========================================================================================
% Running SVD analysis per spectrum:

Fraction=zeros(Ny,Nx,Ns);
Err=zeros(Ny,Nx);

n=0;
Tstart=now;

for j=1:Ny
    for k=1:Nx
        n=n+1;
        S=squeeze(LA(j,k,:));
        fr=UnMix*S;
        Fraction(j,k,:)=fr;                    % Fraction of each component is the output of analysis
        Err(j,k)=std(S-A*fr)^2;
        
        if fix(n/Nr)==n/Nr                     % progress report
            pFinish=100*n/(Nx*Ny);
            tPass= (now-Tstart)*24*60;
            tLeft = tPass*(100-pFinish)/pFinish;
            if ProgRep=='t'
                % progress report by text:
                disp([num2str(pFinish) ' % ||  Time from start: ' num2str(tPass ) ' min. ||  Remaining Time: ' num2str(tLeft)  ' min.  '])
            elseif ProgRep=='g'
                % progress report by graph:
                figure(22); clf;
                C1=A(:,1)*fr(1);
                C2=A(:,2)*fr(2);
                C3=A(:,3)*fr(3);
                C4=A(:,4)*fr(4);
                plot(L,S,'or', L,A*fr,'-b',L,C1,'--c',L,C2,'--c',L,C3,'--c',L,C4,'--c')
                drawnow
                title(['Error (SSE): ' num2str( Err(Ny,Nx)) ' ||  Finished: ' num2str(pFinish) ' % ||  Time from start: ' num2str(tPass ) ' min. ||  Remaining Time: ' num2str(tLeft)  ' min.  ']);
                pause(0.1)
            end
        end
    end
end
end


function [IndAll,Spectra]=ShowSpec2(RGBHE,spec,lambda,n,Fraction)
% RGBHE - RGB Image for the original H&E transmission
% spec - data cube
% labmbda ...
% n - number of nearest neibors for averaging when displaying spectrum.

FRGB=FalseFractionRGB(Fraction,[0,0,3,0; 0,0.7,0,0; 1,0,0,1]);      % false RGB representing SVD analysis results (Fraction)

f=figure(6);
if isempty(f.CurrentAxes)               % if figure exist no need to redraw so magnif remains as before
    subplot(2,2,1); 
    imshow(RGBHE)  
    a=f.CurrentAxes;
    subplot(2,2,2);
    imshow(FRGB)
    b=f.CurrentAxes;
end

a.Position=[0.0    0.48    0.6    0.5];
b.Position=[0.5    0.48    0.6    0.5];

IndAll=[];
Spectra=[];

stop=0;

while stop~=1
    [x,y,button]=ginput(1);
    X=round(x);
    Y=round(y);
    Ind=[X Y];
    if button=='s' || button=='S'       % 's' for stopping
        stop=1; 
    elseif button=='p' || button=='P'   % 'p' may be used for pausing and changing magnification, hit Enter when done 
        pause
        % adjust magnif similar for both views:
        c=f.CurrentAxes;    % read imag limits from last touched subplot
        XL=c.XLim;          
        YL=c.YLim;
        a.XLim=XL;          % apply same limits to both plots
        a.YLim=YL;
        b.XLim=XL;
        b.YLim=YL;
    elseif button=='n' || button=='N'   % 'n' for changing the averaging area
        n=input('Averaging area size? : ');
    else                                % mouse-click - choosing a point a drawing spectrum
        subplot(2,2,3); 
        a3=f.CurrentAxes;
        a3.Position=[0.1300    0.0800    0.3347    0.3];
        if n==0
            Spectrum=squeeze(spec(Y,X,:));
            F1=Fraction(Y,X,1);
            FE=Fraction(Y,X,2);
            FH=Fraction(Y,X,3);
            F4=Fraction(Y,X,4);
        else
            Spectrum=squeeze(mean(spec(Y-n:Y+n,X-n:X+n,:),[1 2]));
            F1=mean(Fraction(Y-n:Y+n,X-n:X+n,1),[1 2]);
            FE=mean(Fraction(Y-n:Y+n,X-n:X+n,2),[1 2]);
            FH=mean(Fraction(Y-n:Y+n,X-n:X+n,3),[1 2]);           
            F4=mean(Fraction(Y-n:Y+n,X-n:X+n,4),[1 2]);           
        end
        plot(lambda,Spectrum)
        axis([400 800 0 1.2])
        title(['X:  ' num2str(X) ' , Y:  ' num2str(Y)])
        IndAll=[IndAll; Ind];
        Spectra=[Spectra; Spectrum];
        
        subplot(2,2,4)
        a4=f.CurrentAxes;
        a4.Position=[0.5703    0.1100    0.3347    0.3412];

        bar(1:4,[F1 FE FH F4])
        axis([0 5 -.2 1.2])
        title(['F1:  ' num2str(F1,'%.3f')   '      , Eosin:  ' num2str(FE,'%.3f') '       ,  Hematoxylin:  ' num2str(FH,'%.3f') '      , F4:  ' num2str(F4,'%.3f') ])
        text(1,1,['n= ' num2str(n)])
  
    end
   
end

end



function Af=SpatialFilter(A,Nf)

% Nf - filter half-size:
x=-Nf:Nf;
y=x';
X=ones(2*Nf+1,1)*x;
Y=y*ones(1,2*Nf+1);
F=exp(-2*(X.^2+Y.^2)/Nf^2);        % Gaussian filter definition
F=F/sum(sum(F));         

Af=conv2(A,F,'same');        % applying filter to Mask
end
function FRGB=FalseFractionRGB(Fraction,W)

% W = weights for RGB false coloring
% e.g. W=[0,0,2,0; 0,0.5,0,0; 1,0,0,1]


FRGB(:,:,1)=W(1,1)*Fraction(:,:,1)+W(1,2)*Fraction(:,:,2)+W(1,3)*Fraction(:,:,3)+W(1,4)*Fraction(:,:,4);      % R
FRGB(:,:,2)=W(2,1)*Fraction(:,:,1)+W(2,2)*Fraction(:,:,2)+W(2,3)*Fraction(:,:,3)+W(2,4)*Fraction(:,:,4);      % G
FRGB(:,:,3)=W(3,1)*Fraction(:,:,1)+W(3,2)*Fraction(:,:,2)+W(3,3)*Fraction(:,:,3)+W(3,4)*Fraction(:,:,4);      % B
end
function XRGB=Spec2RGB(spec,lambda)

R=find(lambda>600);
G=find((lambda>500).* (lambda<600));
B=find(lambda<500);
Dim=size(spec);

XRGB(:,:,1)=sum(spec(:,:,R),3);
Map=squeeze(XRGB(:,:,1));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,1)=XRGB(:,:,1)./prctile(Sig,98);

XRGB(:,:,2)=sum(spec(:,:,G),3);
Map=squeeze(XRGB(:,:,2));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,2)=XRGB(:,:,2)./prctile(Sig,98);

XRGB(:,:,3)=sum(spec(:,:,B),3);
Map=squeeze(XRGB(:,:,3));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,3)=XRGB(:,:,3)./prctile(Sig,98);
end
%}

%-----------------------------------------------------------------

% mask an img on another.
%{
%{
%put kernel mask filtered img on whitened img
mask_kernel = imread('C161VA_x20_X700Y600_20_kmeans_c_201_242_15_3000_on_k13_20_1000.png');
origin_img = imread('C161VA_x20_X700Y600_Whitened_2.png');
masked_img=add_mask_on_img(origin_img, mask_kernel);
imtool(masked_img);
%}
%{
rgb_classic = imread('C161VA_x20_X700Y600_Whitened.png');
mask_kernel = imread('kernel_8_300_on_empty_filtered_img.png');
masked_img=add_mask_on_img(rgb_classic, mask_kernel);
imtool(masked_img);
%}

mask_kernel= imread('C161VA_x20_X700Y600_20_kmeans_c_40_242_200_150000000_on_empty.png');
origin_img = imread('c161_HnE.png');
masked_img=add_mask_on_img(origin_img, mask_kernel);
imtool(masked_img);
%}


%-----------------------------------------------------------------

%this function adds ONE  circle to an img
%{
%{
size_1=(819-520);
size_2=(699-400);
img =load('smthg')
indx = 34;
ct = centers(indx,:);
rd=radii(indx);
create_check_circle_mat(img, ct, rd);
imtoool(final);
%}
function [final] = create_check_circle_mat(img, ctrs, rd)
    %check_mat = zeros(size_1, size_2);

     check_mat = insertShape(img,'FilledCircle', [ctrs(1) ctrs(2) rd], 'Color', 'white','Opacity',1);
     
     final = rgb2gray(check_mat);
     imtool(final);
end
%}

%-----------------------------------------------------------------

%this function adds circles to an img
%{
%{
crop_k=imread('P055A_X2000Y0_R1_Whitened.png');
%temp = imread('SpecToRgb.png');
%crop_k = imcrop(temp,[520 400  (819-520) (699-400)]);
centers_list = centers;
radii_list = radii;

img_with_circles= add_circles_to_img(crop_k, centers_list, radii_list);
%}

%{
function[img_with_circles]= add_circles_to_img(img, centers_list, radii_list)
%this function adds circles to an img

    img_with_circles= img;
    [a,~]=size(centers_list);
    [b,~]=size(radii_list);
    if a~=b
        frpintf('size err: centers and radii list have different sizes. \n');
    end
    for i=1:a
        %fprintf('i= %d, centers : %d, %d , radii: %d \n ', i, centers(i,1), centers(i,2), radii(i));
        img_with_circles = insertShape(img_with_circles,'Circle', [centers_list(i,1) centers_list(i,2) radii_list(i)],'LineWidth',1,'Color','red');
    end
    imtool(img_with_circles)
end
%}
%}
%-----------------------------------------------------------------

%CC tools
%{
CC = bwconncomp(binay);
CC.NumObjects
%CC.PixelIdxList
stats = regionprops('Table',CC,'basic');
centroid_list = stats.Centroid;
ctrd = centroid_list(i,:);

number_of_elt_in_ccs = stats.Area;
BoundingBox = stats.BoundingBox;
fprintf('bb %f \n', BoundingBox(i,:));
[ca, ~] = size(centroid_list);

 %}

%-----------------------------------------------------------------

%KMEANS
%{
%mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/3.rawstack.mat');
%{
mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/C161VA_x20_X700Y600_Whitened_2');

input_spectral_img = mymat.spec(:,:,:);
number_of_clusters = 3;
%1 if want to save img
save_kmean_img=0;
name_of_input_img = 'C161VA_x20_X700Y600_Whitened_2';
%}

%{
noiseless = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/c161va_without_noise.mat');
input_spectral_img = noiseless.FImage;
number_of_clusters = 20;
save_kmean_img=0;
name_of_input_img = 'noiseless';
%}

%{

noiseless = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/P055A_X2000Y0_R1.mat');
input_spectral_img = noiseless.spec;
number_of_clusters = 20;
save_kmean_img=0;
name_of_input_img = 'p055';
%}


%{
median_smooth = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/smoothed_median_C161va_x20_x700y600.mat');
mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/C161VA_x20_X700Y600.mat');
wl_mat =mymat.spec;
input_spectral_img =median_smooth.new_mat;
%wl_mat(1:28);
number_of_clusters = 2;
save_kmean_img=0;
name_of_input_img = 'smoothed_median_c161_x20';
%}

kmean_mat_output=kmeans_algo(input_spectral_img, number_of_clusters, save_kmean_img, name_of_input_img);
imtool(kmean_mat_output);

%}

%-----------------------------------------------------------------

%filtered_selection
%{
%parameters
% SELECT working_img : which kmeans to work on . 
%{
%temp = imread('10_kmeans_of_mat_P055A_X2000Y0_R1.png');
%crop_k =mat2gray(temp);

%first img loads
temp = imread('10K_means_out.png');
cropk= imcrop(temp,[520 400  (819-520) (699-400)]);
%imtool(cropk);
%}

%reference cell
wkn_img = imread('C161VA_x20_X700Y600_20_kmeans.png');
%crop_k=imcrop(wkn_img,[520 400  (819-520) (699-400)]);
name = 'C161VA_x20_X700Y600_20_kmeans';

% select  clusters to filter 
%{
%first img clusters values
first_kernel = [0,28,142,198,170];
first_cytoplasm = [227,255];
first_membrane = [55,113,57];
%}

%reference cell clusters values

%C161VA_x20_X700Y600_20_kmeans
k_13 =[13];
reference_kernel = [13,121, 242, 40];
reference_cyto = [201,161, 148];
c_201_242 = [201,242];
c_40 = [40];
c_134 = [134];
c_40_242 = [40,242];
c_148 = [148];
temp  =[ 40, 242];
cancerous_zone = [40,134 ,201 , 148,242];
canc_cell =[.16];
cell_area = c_40_242;
cell_area_name = 'c_40_242';


% select min - max allowed size
Max_allowed_cc_size = 15;
Min_allowed_cc_size = 0;


%select background
%{
%kernel_out_8_300 = imread('kernel_out_8_300.png');
%support_name= 'kernel_8_300';
%support=mat2gray(kernel_out_8_300);

%m_k_8_300 = imread('kernel_8_300_on_empty_filtered_img.png');
%support_name= 'm_k_8_300';
%support=mat2gray(m_k_8_300);

%membrane_kernel_8_100 = imread('kernel_out_8_100.png');
%support =mat2gray(membrane_kernel_8_100);
%support_name= 'membrane_kernel_8_100';


[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;


background = mat2gray(imread('C161VA_x20_X700Y600_20_kmeans_k_13_20_1000_on_empty.png'));
bg_name= 'k13_20_1000';

%}
[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;



bg_reduction = 0;
save_filtered_img = 0;


final_name = append(name, '_', cell_area_name,'_',int2str(Min_allowed_cc_size),'_',int2str(Max_allowed_cc_size),'_on_',bg_name);

f_i = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction);
imtool(f_i);

%{
kernel = -1*ones(3);
kernel(2,2) = 17;
enhancedImage = imfilter(f_i, kernel);
enhancedImage = imfill(enhancedImage,'holes');
imtool(enhancedImage);
imwrite(enhancedImage,'C161VA_x20_X700Y600_cancerous_zone.png');
%}

%}

%-----------------------------------------------------------------

%SPECTORGB _ LIOR
%{
%{
mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/3.rawstack.mat');
mat = mymat.spec(:,:,:);
lambda = mymat.lambda; 
rgbImg= spectoRGB(mat,lambda, 0,'whatever');
imtool(rgbImg);
%}
%{
mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/C161VA_x20_X700Y600');
save_SpecToRgb=0;
name= 'C161VA_x20_X700Y600';
lambda = mymat.lambda;
mat = mymat.spec(:,:,:);
rgbImg= spectoRGB(mat, lambda, save_SpecToRgb,name);
imtool(rgbImg);
%}
%{
mymat= load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/EmptyRef.mat');
save_SpecToRgb=0;
name= 'empty_ref';
lambda = mymat.lambda;
mat = mymat.spec(:,:,:);
rgbImg= spectoRGB(mat,lambda,save_SpecToRgb,name);
imtool(rgbImg);
%}

%{
noiseless = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/c161va_without_noise.mat');
input_spectral_img = noiseless.FImage;
%imshow(input_spectral_img);
save_SpecToRgb=0;
name= 'noiseless_c161va_test_1';
empty_ref= load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/EmptyRef.mat');
lambda = empty_ref.lambda;
mat = input_spectral_img;
rgbImg= spectoRGB(mat,lambda,save_SpecToRgb,name);
imtool(rgbImg);
%}

%{
%with smoothing wl

smoothed_median = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/smoothed_median_uint8_C161va_x20_x700y600.mat');
%input_spectral_img =  cast(smoothed_median.new_mat, 'uint8'); 
%imshow(input_spectral_img);
save_SpecToRgb=0;
name= 'smoothed_median_c161va';
empty_ref= load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/EmptyRef.mat');
mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/C161VA_x20_X700Y600.mat');
lambda = mymat.lambda(1:28);
mat = double(smoothed_median.new_mat);
rgbImg= spectoRGB(mat,lambda,save_SpecToRgb,name);
imtool(rgbImg);
%}
%}


%-----------------------------------------------------------------

%whitening rgb img
%{
%{
rgbImage = imread('SpecToRgb.png');

top_left_corner_x_area1 =946;
top_left_corner_y_area1= 609;
bottom_right_corner_x_area1 =982;
bottom_right_corner_y_area1 =639;

top_left_corner_x_area2 =597;
top_left_corner_y_area2= 736;
bottom_right_corner_x_area2 =621;
bottom_right_corner_y_area2 =754;

save_new_white=0;
name = 'P055A_X2000Y0_R1';
%}

%whitening P055A_X2000Y0_R1 img
rgbImage = imread('C161VA_x20_X700Y600_SpecToRgb.png');
top_left_corner_x_area1 =1110;
top_left_corner_y_area1= 292;
bottom_right_corner_x_area1 =1210;
bottom_right_corner_y_area1 =356;

top_left_corner_x_area2 =1962;
top_left_corner_y_area2= 276;
bottom_right_corner_x_area2 =2060;
bottom_right_corner_y_area2 =383;

save_new_white=0;
name = 'C161VA_x20_X700Y600';


imtool(rgbImage);
new_white = whitening_rgb (rgbImage,top_left_corner_x_area1,top_left_corner_y_area1,bottom_right_corner_x_area1,bottom_right_corner_y_area1,top_left_corner_x_area2,top_left_corner_y_area2,bottom_right_corner_x_area2,bottom_right_corner_y_area2,save_new_white,name);
imtool(new_white);
%}

%-----------------------------------------------------------------

%smothing and cutting wl from spec
%{


%some general filters

horiz_filter = [0 0 0; 1 1 1; 0 0 0];
vert_filter = [0 1 0; 0 1 0; 0 1 0];
ones_filt = [ 1 1 1; 1 1 1; 1 1 1];
mean_filter = ones(3,3)/9;
%meadian_filter = medfilt2(img);
gaussian = fspecial('gaussian');
laplacian = fspecial('laplacian');

%remind step3 remove noise
%{
bw = bwareaopen(bw,30);
%}


%remind Step 4: Find the Boundaries
%{
[B,L] = bwboundaries(bw,'noholes');

imshow(label2rgb(L,@jet,[.5 .5 .5]))
hold on
for k = 1:length(B)
  boundary = B{k};
  plot(boundary(:,2),boundary(:,1),'w','LineWidth',2)
end
%}
mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/C161VA_x20_X700Y600.mat');
wl_mat =mymat.spec;
lambda_mat = mymat.lambda;

chosen_wl=15;
grayscale_out = read_unique_wl(wl_mat , lambda_mat, chosen_wl);
%imtool(grayscale_out);

%smoothing
smoothed = imfilter(grayscale_out,mean_filter, 'same');
%or conv2 or filter2 for conv
imtool(uint8(255 * mat2gray(smoothed)));



%  same smoothing to all WL , and  save the new spectrum
%{
[a,b] = size( grayscale_out);
smoothed= medfilt2(grayscale_out);
new_mat = smoothed;

for i =2:28
chosen_wl=i;
grayscale_out = read_unique_wl(wl_mat , lambda_mat, chosen_wl);
smoothed= medfilt2(grayscale_out);
new_mat = cat(3,new_mat, smoothed);
end
size(new_mat)
%new_mat  = cast(new_mat, 'double');
name = "smoothed_median_uint8_C161va_x20_x700y600.mat";
save(name,'new_mat');
%}

%}

%-----------------------------------------------------------------

%edge detect
%{
median_smooth_rgb = imread('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/smoothed_median_c161va_SpecToRgb.png');
I=rgb2gray(median_smooth_rgb);
edges = edge(I,'sobel',0.015);
%imwrite(edges, 'c161va_edge_sobel_.015.png');
imtool(edges);
%}

%-----------------------------------------------------------------

%beside
%{
%unet from file
%https://github.com/clamesc/Machine-Learning-in-Medical-Imaging--U-Net

%draw circles:
%viscircles(centers,radii)
%imtool(first_circle_out)


%website for detecting round + roundness
%{

%step1 : load image
RGB = imread('pillsetc.png');
%imshow(RGB)


%Step 2: Threshold the Image
I = rgb2gray(RGB);
bw = imbinarize(I);
%imshow(bw)

%imshow(bw)

%fill gap 
se = strel('disk',2);
bw = imclose(bw,se); 
bw = imfill(bw,'holes');


%Step 4: Find the Boundaries
[B,L] = bwboundaries(bw,'noholes');

imshow(label2rgb(L,@jet,[.5 .5 .5]))
hold on
for k = 1:length(B)
  boundary = B{k};
  plot(boundary(:,2),boundary(:,1),'w','LineWidth',2)
end

%Step 5: Determine which Objects are Round
stats = regionprops(L,'Area','Centroid');

threshold = 0.94;

% loop over the boundaries
for k = 1:length(B)

  % obtain (X,Y) boundary coordinates corresponding to label 'k'
  boundary = B{k};

  % compute a simple estimate of the object's perimeter
  delta_sq = diff(boundary).^2;    
  perimeter = sum(sqrt(sum(delta_sq,2)));
  
  % obtain the area calculation corresponding to label 'k'
  area = stats(k).Area;
  
  % compute the roundness metric
  metric = 4*pi*area/perimeter^2;
  
  % display the results
  metric_string = sprintf('%2.2f',metric);

  % mark objects above the threshold with a black circle
  if metric > threshold
    centroid = stats(k).Centroid;
    plot(centroid(1),centroid(2),'ko');
  end
  
  text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','y',...
       'FontSize',14,'FontWeight','bold')
  
end

title(['Metrics Closer to 1 Indicate that ',...
       'the Object is Approximately Round'])

%end of steps, trying to find circle using imfindcircles
%[centers, radii, metric] = imfindcircles(rgb,[10 30],'ObjectPolarity','dark', 'Sensitivity', .97);
%centers
%figure
%imshow(rgb)
%h = viscircles(centers,radii);
%}



%lior and franch tag CC for course ibud tmuna
%{
function [newImg] = tagConnectedComponent(img)
%TAGCONNECTEDCOMPONENT Lior shimon 341348498 Franck schwartz 329863237 . this function transform the input img into a binary img,then use mat_tranforn function to get all the potential connected componnent, iterate over it and finally use conv_mat function to output the final binary connected component matrix. 
count = 0;
i_binary = imread(img);
%meanIntensity = mean(i_binary(:));
%i_binary = i_binary > meanIntensity;

%i_binary = [1 0 1 0 0 0 1; 1 0 1 0 1 0 1; 1 0 1 0 1 0 1; 1 1 1 1 1 0 1; 0 0 0 0 0 0 1; 0 0 1 1 1 1 1];
%[r,c] = size(i_binary);
[mat_Tag, save_label]= mat_transform(i_binary);
fprintf('mat_transform function terminated \n');
mat_neigh = save_label;
keep = 1;
[r,c] = size(save_label);
temp = zeros(r,c);
log_size = log(max(r,c));
while keep == 1 || count <= log_size
    mat_neigh = mat_neigh * save_label;
    mat_neigh = double(mat_neigh > 0);
    count = count + 1;
    if mat_neigh == temp
        keep = 0;
    end
    temp = mat_neigh;
end
CV = conv_mat(mat_neigh,c);
fprintf('conv_mat function terminated \n');

newImg = tag_comp(mat_Tag, CV);
fprintf('tag_comp function terminated \n');

%disp(newImg);
%imshow(newImg);
end

function [tag_mat, save_lab] = mat_transform(img)
%MAT_TRANSFORM this function run over all the pixels, and give a specificlabel to each potential connected component . It returns the "tagged matrix" and save_lab will be used to create the convertor vector later
 imgTemp = img;
[row,col] = size(img);
left = zeros(row+1,1);
up = zeros(1,col);
imgTemp = [up; imgTemp];
imgTemp = [left imgTemp];
labels = 1;
diag_label = ones(1,col*3);
save_label = diag(diag_label);
for r = 2:row + 1
        for c = 2:col + 1
            if imgTemp(r,c) ~= 0
                if imgTemp(r-1,c) ~= 0 
                    if imgTemp(r,c-1) ~= 0 
                        save_label(imgTemp(r-1,c),imgTemp(r,c-1)) = 1;
                        save_label(imgTemp(r,c-1), imgTemp(r-1,c)) = 1;
                    end
                    imgTemp(r,c) = imgTemp(r-1,c);
                elseif imgTemp(r,c-1) ~= 0
                    imgTemp(r,c) = imgTemp(r,c-1);
                else
                    imgTemp(r,c) = labels;
                    labels = labels + 1;
                end
            end
        end
end
lab = labels -1;
save_lab = save_label(1:lab, 1:lab);
%disp(save_lab);
tag_mat = imgTemp(2:row+1, 2:col+1);
%disp(tag_mat);
end

function [vec] = conv_mat(mat_neigh, n)
%CONV_MAT this function crete the converter vector from the labelled Marix
%in the way we learnt in the tirgul
CV = mat_neigh(1,:);
label = 2;

for i=2:n
    if CV(i) == 0
        CV = CV + mat_neigh(i,:) * label;
        %disp(CV);
        label = label + 1;
    end
end
vec = CV;
%disp(vec);
end

function [mat_lab] = tag_comp(mat_tag, CV)
%TAG_COMP this function take the converter vector and the  mat_tag matrix and
%create the sequenced connected component, i.e the final expected matrix.
imgTmp = mat_tag;
[row,col] = size(mat_tag);
for r = 1:row
        for c = 1:col
            if imgTmp(r,c) ~= 0
                imgTmp(r,c) = CV(imgTmp(r,c));
            end
        end
end
mat_lab = imgTmp;
%disp(mat_lab);
end
%}
%}

%-----------------------------------------------------------------
%retrove specific function
%{

%{
tissue_name ='x40_r1_r2_spec';
struct_path = '/Users/lior/Desktop/Image & Analysis /' +  string(tissue_name) + '.mat';
data = load(struct_path);
%RGBHE = data.rgb;
wl_mat = data.spec;
%FImage = immg = wl_mat
lambda_mat = data.lambda;
%lambda = lambda_mat
%}
%{
function [Normalized_H, Normalized_E] = get_HnE(wl_mat,lambda_mat)
[Fraction,Err]=SVDAnalysisHnE(wl_mat,lambda_mat,'g',1);
E= Fraction(:,:,2);
Normalized_E = uint8(255 * mat2gray(E));
%imtool(Normalized_E);
H= Fraction(:,:,3);
Normalized_H = uint8(255 * mat2gray(H));
%imtool(Normalized_H);




W1=[0,0,2,0; 0,0.5,0,0; 1,0,0,1];
W2=[0,0,3,0; 0,0.7,0,0; 1,0,0,1];


FRGB=FalseFractionRGB(Fraction,W1);

imtool(uint8(255 * mat2gray(FRGB(:,:,1))));
imtool(uint8(255 * mat2gray(FRGB(:,:,2))));

imtool(FRGB);


tissue_name ='x40_r1_r2_spec';
struct_path = '/Users/lior/Desktop/Image & Analysis /' +  string(tissue_name) + '.mat';
data = load(struct_path);

rgb = data.rgb;
lambda_mat = data.lambda;
wl_mat = data.spec;

FImage = wl_mat;
lambda = lambda_mat;
[Fraction,Err]=SVDAnalysisHnE(FImage,lambda,'t',0);


FRGB=FalseFractionRGB(Fraction,W1);      % false RGB representing SVD analysis results (Fraction)
imtool(FRGB);


E= Fraction(:,:,2);
nE = uint8(255 * mat2gray(RGBHE));

H= Fraction(:,:,3);


%imtool(E);






%imtool(nE);
%imwrite(nE, 'current_work_e.png');

%imtool(H);

nH = uint8(255 * mat2gray(RGBHE));
imtool(nH);
q=1;
%imwrite(nH, 'current_work_h.png');
%}

%}

%-----------------------------------------------------------------
