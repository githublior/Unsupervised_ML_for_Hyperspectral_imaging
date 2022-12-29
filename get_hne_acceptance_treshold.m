function value = get_hne_acceptance_treshold( mat, qtle);
%quantile : [0,1]

mat = reshape(mat.',1,[]);
value = quantile(mat, qtle);
end

