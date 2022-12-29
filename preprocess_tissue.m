function [data, hne, whitened_img ] = preprocess_tissue(tissue_struct, whitening)

%default cut in spectrum dimension:
bottom_wl_treshold = 4;
top_wl_treshold = 30;


rgb = tissue_struct.RGB;
wl_mat = tissue_struct.spec;
lambda_mat = tissue_struct.lambda;


%-----------------------------------------------------------------

%##############  Cutting and smothing  wl from spec
 %compute  cutting and smoothing
smoothed_and_wl_shaped = medfilt2(read_unique_wl(wl_mat , lambda_mat, bottom_wl_treshold));
for i = bottom_wl_treshold+1 :top_wl_treshold
single_wl = read_unique_wl(wl_mat , lambda_mat, i);
smoothed_wl= medfilt2(single_wl);
smoothed_and_wl_shaped = cat(3,smoothed_and_wl_shaped, smoothed_wl);
end
wl_mat = smoothed_and_wl_shaped;
lambda_mat = lambda_mat(bottom_wl_treshold:top_wl_treshold);


%VISU
%{
cutted_and_smoothed_spectral_img =  spectoRGB(smoothed_and_wl_shaped,lambda_mat, 0,'name');
cutted_and_smoothed_spectral_img_formatted =  uint8(255 * mat2gray(cutted_and_smoothed_spectral_img));
imshow(cutted_and_smoothed_spectral_img_formatted)
%}
%SAVING: 
%{
%inutile%smoothed_and_wl_shaped  = cast(new_mat, 'double');
name = string(tissue_name )+ "_smoothed_and_wl_shaped.mat" ;
save(name,'wl_mat', 'lambda_mat', 'rgb');
%}



%OR load cutting and smoothing
%{
data = load(fullfile(pwd,tissue_name + ".mat" ));
rgb = data.RGB;
wl_mat = data.spec;
lambda_mat = data.lambda;
%}



data.RGB = rgb;
data.wl_mat = wl_mat;
data.lambda_mat = lambda_mat;

%-----------------------------------------------------------------


%Compute HnE
[H, E] = get_HnE(wl_mat,lambda_mat);
hne.H = H;
hne.E = E;



%-----------------------------------------------------------------

if whitening
    fprintf('select top left corner white area 1');
    imshow(rgb);
    [xi, yi, ~] = ginput(1); 
    top_left_corner_x_area1 = uint8(xi);
    top_left_corner_y_area1 = uint8(yi);
    fprintf('select bottom right corner white area 1');
    imshow(rgb);
    [xi, yi, ~] = ginput(1); 
    bottom_right_corner_x_area1 = uint8(xi);
    bottom_right_corner_y_area1 =  uint8(yi);
    
    
    fprintf('select top left corner white area 2');
    imshow(rgb);
    [xi, yi, ~] = ginput(1); 
    top_left_corner_x_area2 =  uint8(xi);
    top_left_corner_y_area2 = uint8(yi);
    fprintf('select bottom right corner white area 2');
    imshow(rgb);
    [xi, yi, ~] = ginput(1); 
    bottom_right_corner_x_area2 =  uint8(xi);
    bottom_right_corner_y_area2 =  uint8(yi);
    
    save_new_white = 0;
    name = ' ';
    
    [whitened_img]= whitening_rgb (rgb,top_left_corner_x_area1,top_left_corner_y_area1,bottom_right_corner_x_area1,bottom_right_corner_y_area1,top_left_corner_x_area2,top_left_corner_y_area2,bottom_right_corner_x_area2,bottom_right_corner_y_area2,save_new_white,name);
else 
    whitened_img = rgb;
end



end