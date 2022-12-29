function distance_map_list = postprocess_on_cluster_selection(relevant_cluster, nbr_of_cluster, D, W_kmean_mat_output)



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


end