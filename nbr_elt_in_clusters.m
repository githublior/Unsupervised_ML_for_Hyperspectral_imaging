function [ L_map ] = nbr_elt_in_clusters(a,b, L) ;
% count # pixels in each CC.
    L_map = containers.Map('KeyType','double','ValueType','double');

    for i = 1:a
        for j= 1:b
            if L(i,j) ~= 0
                if isKey(L_map,L(i,j))
                     L_map(L(i,j))=  L_map(L(i,j)) +1;
                else
                    L_map(L(i,j))= 1;
                    %fprintf('new value detected: %d \n', L(t1,t2))
                end
            end
        end
    end
