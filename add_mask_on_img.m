function [masked_img]=add_mask_on_img(origin_img, mask_kernel);
%mask is 0-1 matrix.
%{
last_1 = origin_img(:, :, 1) .* (1-mask_kernel);
last_2 = origin_img(:, :, 2) .* (1-mask_kernel);
last_3 = origin_img(:, :, 3) .*  (1-mask_kernel);

masked_img = uint8(zeros(size(origin_img)));
%imtool(last_1);
masked_img(:,:,1)= uint8(last_1);
masked_img(:,:,2)= uint8(last_2);
masked_img(:,:,3)= uint8(last_3);
%}
masked_img(:,:,1)=origin_img(:,:,1).*(1-mask_kernel);
masked_img(:,:,2)= origin_img(:,:,2).*(1-mask_kernel);
masked_img(:,:,3)= origin_img(:,:,3).*(1-mask_kernel);