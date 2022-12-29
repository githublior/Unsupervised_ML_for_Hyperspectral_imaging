%{
parametre classique pour utiliser cette function
tissue_name ='x40_r1_r2_spec';
struct_path = '/Users/lior/Desktop/Image & Analysis /' +  string(tissue_name) + '.mat';
data = load(struct_path);
%RGBHE = data.rgb;
wl_mat = data.spec;
%FImage = immg = wl_mat
lambda_mat = data.lambda;
%lambda = lambda_mat
%}

function [Normalized_H, Normalized_E] = get_HnE(wl_mat,lambda_mat)
[Fraction,~]=SVDAnalysisHnE(wl_mat,lambda_mat,'n',1);
E= Fraction(:,:,2);
Normalized_E = uint8(255 * mat2gray(E));
H= Fraction(:,:,3);
Normalized_H = uint8(255 * mat2gray(H));
end
%imtool(Normalized_H);
%imtool(Normalized_E);