function [origin_img]=red_mask_on_img(origin_img, mask_kernel);
%mask is 0-1 matrix.


[a,b,~] = size(mask_kernel);
for i=1:a
    for j=1:b
        %fprintf('%d', mask_kernel(i, j));
       
        if mask_kernel(i, j) >0
            origin_img(i, j, 1)= 255;
            origin_img(i, j, 2)= 0;
            origin_img(i, j, 3)= 0;
        end
        
    end
end
