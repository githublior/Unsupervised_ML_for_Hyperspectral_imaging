%SPECTORGB _ LIOR

%{
%{
mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/3.rawstack.mat');
mat = mymat.spec(:,:,:);
lambda = mymat.lambda; 
rgbImg= spectoRGB(mat,lambda, 0,'no');
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


function rgbImg= spectoRGB(mymat,lambda, save_SpecToRgb,name)
%{
- mymat : mat of spectral img.
 lambda- wavelength list associated to mymat values.
-save_SpecToRgb: boolean. 1 if you want to save the output, 0 else
-name : string. name when saving output
%}

[a,b,c] = size(mymat);
Y = reshape(mymat, [a*b, c]);
sRGB = spectrumRGB(lambda); 
size(sRGB);
X = reshape(sRGB, [1*c,3]);
size(X)
size(Y)
z =Y*X;
Z = reshape(z,[a,b, 3]);

rgbImg =  uint8(255 * mat2gray(Z));

if save_SpecToRgb==1
    nm =append(name,'_SpecToRgb', '.png' );
    imwrite(rgbImg, nm);
end
%imshow(rgbImg);


fprintf('spectoRGB function terminated. \n');
end


