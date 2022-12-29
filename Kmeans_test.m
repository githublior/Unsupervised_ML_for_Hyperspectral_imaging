
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
number_of_clusters = 20;
save_kmean_img=1;
name_of_input_img = 'smoothed_median_c161_x20';
%}

kmean_mat_output=kmeans_algo(input_spectral_img, number_of_clusters, save_kmean_img, name_of_input_img);
imtool(kmean_mat_output);

%}
function kmean_mat_output=kmeans_algo(input_spectral_img, number_of_clusters, save_kmean_img, name_of_input_img)
%{
-input_spectral_img: input matrix
- number_of_clusters : K
-save_kmean_img : boolean. 1 to save the output, 0 else
-name_of_input_img name to use when saving output
%}

    K=number_of_clusters;
    rng(1); % For reproducibility
    [a,b,c]=size(input_spectral_img);

    %reshape input_img for applying kmeans on a 2d matrix.
    X = reshape(input_spectral_img, [a*b, c]);

    % KMeans algo
    [idx2Region,~] = kmeans(X,K);

    %put the clusterization on a new img
    A = zeros(a*b,1);
    for i = 1:a*b
        %{
        if (i > 1087*500) && (i < 1087*600)
            A(i) = K-1;
        else
            A(i) = idx2Region(i,1)-1;
        end
        %}
        A(i) = idx2Region(i,1)-1;
    end

    A = reshape(A, [a,b]);
    kmean_mat_output = mat2gray(A,[0,K-1]);

    if save_kmean_img==1
        name =append(name_of_input_img,'_',int2str(K), '_kmeans.png' );
        imwrite(kmean_mat_output, name);
    end

    %figure
    %legend('Region 1','Region 2','Region 3','Data','Location','SouthEast');
   
    
    fprintf('kmeans_algo function terminated. \n');
end

 


 