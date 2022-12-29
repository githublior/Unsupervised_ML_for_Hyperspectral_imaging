
%whitening rgb img
%{
%{
rgbImage = imread('SpecToRgb.png');

top_left_corner_x_area1 =946;
top_left_corner_y_area1= 609;
bottom_right_corner_x_area1 =982;
bottom_right_corner_y_area1 =639;

top_left_corner_x_area2 =597;
top_left_corner_y_area2= 736;
bottom_right_corner_x_area2 =621;
bottom_right_corner_y_area2 =754;

save_new_white=0;
name = 'P055A_X2000Y0_R1';
%}

%whitening P055A_X2000Y0_R1 img
rgbImage = imread('C161VA_x20_X700Y600_SpecToRgb.png');
top_left_corner_x_area1 =1110;
top_left_corner_y_area1= 292;
bottom_right_corner_x_area1 =1210;
bottom_right_corner_y_area1 =356;

top_left_corner_x_area2 =1962;
top_left_corner_y_area2= 276;
bottom_right_corner_x_area2 =2060;
bottom_right_corner_y_area2 =383;

save_new_white=0;
name = 'C161VA_x20_X700Y600';


imtool(rgbImage);
new_white = whitening_rgb (rgbImage,top_left_corner_x_area1,top_left_corner_y_area1,bottom_right_corner_x_area1,bottom_right_corner_y_area1,top_left_corner_x_area2,top_left_corner_y_area2,bottom_right_corner_x_area2,bottom_right_corner_y_area2,save_new_white,name);
imtool(new_white);
%}


function [new_white]= whitening_rgb (rgbImage,top_left_corner_x_area1,top_left_corner_y_area1,bottom_right_corner_x_area1,bottom_right_corner_y_area1,top_left_corner_x_area2,top_left_corner_y_area2,bottom_right_corner_x_area2,bottom_right_corner_y_area2,save_new_white,name)
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
end



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


%make the means of each color channel the same
%{
grayImage = rgb2gray(rgbImage); % Convert to gray so we can get the mean luminance.
% Extract the individual red, green, and blue color channels.
redChannel = rgbImage(:, :, 1);
greenChannel = rgbImage(:, :, 2);
blueChannel = rgbImage(:, :, 3);
meanR = mean2(redChannel);
meanG = mean2(greenChannel);
meanB = mean2(blueChannel);
meanGray = mean2(grayImage);
% Make all channels have the same mean
redChannel = uint8(double(redChannel) * meanGray / meanR);
greenChannel = uint8(double(greenChannel) * meanGray / meanG);
blueChannel = uint8(double(blueChannel) * meanGray / meanB);
% Recombine separate color channels into a single, true color RGB image.
rgbImage = cat(3, redChannel, greenChannel, blueChannel);
%imshow(rgbImage);
%}
%{
%handed white

Y = reshape(rgbImage, [1087*1479, 3]);
for i = 2:1087*1479
    Y(i,1) =  Y(i,1)     ;
    Y(i,2) =  Y(i,2)   +80;%111
    Y(i,3) =  Y(i,3)  + 0     +160; %219
end

rgbImage = reshape(Y,[1087,1479, 3]);


%save and show result
imwrite(rgbImage, 'whiten.png');
%temp = imread('whiten.png');
imtool('whiten.png');

%}
