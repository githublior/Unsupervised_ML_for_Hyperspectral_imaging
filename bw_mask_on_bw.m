function [origin_img]=bw_mask_on_bw(origin_img, mask_kernel);
%mask is 0-1 matrix.


[a,b,~] = size(mask_kernel);
for i=1:a
    for j=1:b
        %fprintf('%d', mask_kernel(i, j));
       
        if mask_kernel(i, j) >0
            origin_img(i, j)= 1;
        end
        
    end
end
