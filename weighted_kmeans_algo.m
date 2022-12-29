function [kmean_mat_output, D] = weighted_kmeans_algo(input_spectral_img, number_of_clusters, save_kmean_img, name_of_input_img,file_path);
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
    [idx2Region,~,~,D]= kmeans(X,K);

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
        name =append(name_of_input_img,'_',int2str(K), '_kmeans_and_D.mat' );
        save(fullfile(file_path,name) , 'kmean_mat_output', 'D');
    end

    %figure
    %legend('Region 1','Region 2','Region 3','Data','Location','SouthEast');
    
    fprintf('distance kmeans_algo function terminated. \n');


 