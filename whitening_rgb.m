function [new_white]= whitening_rgb (rgbImage,top_left_corner_x_area1,top_left_corner_y_area1,bottom_right_corner_x_area1,bottom_right_corner_y_area1,top_left_corner_x_area2,top_left_corner_y_area2,bottom_right_corner_x_area2,bottom_right_corner_y_area2,save_new_white,name);
    %{
    This function correct the colors of the img by selecting 2 rectangles
    defined as white area, and facorize relatively each pixels of the img
    according to the actual tent.
    %}

    %get the white area
    focus_area_1 = imcrop(rgbImage, [top_left_corner_x_area1 top_left_corner_y_area1 (bottom_right_corner_x_area1-top_left_corner_x_area1) (bottom_right_corner_y_area1-top_left_corner_y_area1)]);
    %[sfa,sfb,trois] = size(focus_area);
    focus_area_2= imcrop(rgbImage, [top_left_corner_x_area2 top_left_corner_y_area2 (bottom_right_corner_x_area2-top_left_corner_x_area2) (bottom_right_corner_y_area2-top_left_corner_y_area2)]);
    %imtool(focus_area);

    %get mean value for each channel for first area
    redChannel_1 = focus_area_1(:, :, 1);
    greenChannel_1 = focus_area_1(:, :, 2);
    blueChannel_1 = focus_area_1(:, :, 3);
    meanR_1 = mean2(redChannel_1);
    meanG_1 = mean2(greenChannel_1);
    meanB_1 = mean2(blueChannel_1);

    %get mean value for each channel for second area
    redChannel_2 = focus_area_2(:, :, 1);
    greenChannel_2 = focus_area_2(:, :, 2);
    blueChannel_2 = focus_area_2(:, :, 3);
    meanR_2 = mean2(redChannel_2);
    meanG_2 = mean2(greenChannel_2);
    meanB_2 = mean2(blueChannel_2);

    %get average for each channel
    meanR = (meanR_1 + meanR_2)/2;
    meanG = (meanG_1 + meanG_2)/2;
    meanB = (meanB_1 + meanB_2)/2;

    new_white = rgbImage;
    %[sfa,sfb,trois]= size(new_white);
    [a,b,~] = size(rgbImage);

    for i=1:a
        for j=1:b
            new_white(i,j,1)= rgbImage(i,j,1)*(255/meanR);
            new_white(i,j,2)= rgbImage(i,j,2)*(255/meanG);
            new_white(i,j,3)= rgbImage(i,j,3)*(255/meanB);
        end
    end


    if save_new_white==1
        nm =append(name,'_Whitened.png' );
        imwrite(new_white, nm);
    end    
    fprintf('whitening_rgb function terminated \n');




%{

Whitener_rgb(temp);
function [whtImg]= Whitener_rgb(mat)
Y = reshape(mat, [1087*1479, 3]);
for i = 2:1087*1479
    Y(i,1) =  Y(i,1)  + 39;
    Y(i,2) =  Y(i,2)  + 50;%111
    Y(i,3) =  Y(i,3)  + 119; %219
end

Z = reshape(Y,[1087,1479, 3]);
whtImg =  uint8(Z);
imshow(whtImg);
end
%}

% blue band from j=225:240 all over i .
