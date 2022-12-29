function relevant_cluster = get_user_cluster_selection(W_kmean_mat_output);
fprintf('select cluster value relevant to cells:\n');
imshow(W_kmean_mat_output);
confirm_cluster_selection  = 'a';
while confirm_cluster_selection ~= 'f'
    pixel_values = impixel;
    confirm_cluster_selection = input(' press "f" on terminal to finish the cluster selection. \n press anything else to restart cluster selction....\n','s');
end
    relevant_cluster =rmmissing(unique(pixel_values(:,1))).';
end
  