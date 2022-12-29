
%{ 
Find CIRCLE in Image sources 

https://blogs.mathworks.com/steve/2012/09/04/detecting-circular-objects-in-images/?s_tid=answers_rc2-2_p5_MLT

%https://fr.mathworks.com/help/images/detecting-a-cell-using-image-segmentation.html
step 1 to 8

https://fr.mathworks.com/matlabcentral/answers/305181-detection-of-elliptical-rings-from-microscopy-images
After this, you can play with other morphological operations to remove disturbances and sharpen the edges. The function bwmorph with the skeletonize operation (documentation here) might help in thinning out the lines.

%https://fr.mathworks.com/help/images/identifying-round-objects.html

%}

%https://fr.mathworks.com/help/images/identifying-round-objects.html
temp = imread('cytoplasm_out_50_3000from_kernel_8_100_.png');
se = strel('disk',2);
bw = imclose(temp,se);
bw = imfill(bw,'holes');

[ctrs, rs] = imfindcircles(rgb,[17 28],'ObjectPolarity','dark');
imshow(temp);
h = viscircles(ctrs,rs);

