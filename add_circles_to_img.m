function[img_with_circles]= add_circles_to_img(img, centers_list, radii_list);
%this function adds circles to an img

    img_with_circles= img;
    [a,~]=size(centers_list);
    [b,~]=size(radii_list);
    if a~=b
        frpintf('size err: centers and radii list have different sizes. \n');
    end
    for i=1:a
        %fprintf('i= %d, centers : %d, %d , radii: %d \n ', i, centers(i,1), centers(i,2), radii(i));
        img_with_circles = insertShape(img_with_circles,'Circle', [centers_list(i,1) centers_list(i,2) radii_list(i)],'LineWidth',1,'Color','red');
    end
