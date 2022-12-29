%Schwartz Franck
%329863237
%Shimon Lior
%341348498
%this function transform the input img into a binary img,
%then use mat_tranforn function to get all the potential connected componnent, 
%iterate over it and finally use conv_mat function to output the final
%binary connected component matrix.



tagConnectedComponent('bigMozart.tiff');
function [newImg] = tagConnectedComponent(img)
%TAGCONNECTEDCOMPONENT Lior shimon 341348498 Franck schwartz 329863237 . this function transform the input img into a binary img,then use mat_tranforn function to get all the potential connected componnent, iterate over it and finally use conv_mat function to output the final binary connected component matrix. 
count = 0;
i_binary = imread(img);
meanIntensity = mean(i_binary(:));
i_binary = i_binary > meanIntensity;

%i_binary = [1 0 1 0 0 0 1; 1 0 1 0 1 0 1; 1 0 1 0 1 0 1; 1 1 1 1 1 0 1; 0 0 0 0 0 0 1; 0 0 1 1 1 1 1];
%[r,c] = size(i_binary);
[mat_Tag, save_label]= mat_transform(i_binary);
fprintf('mat_transform function terminated \n');
mat_neigh = save_label;
keep = 1;
[r,c] = size(save_label);
temp = zeros(r,c);
log_size = log(max(r,c));
while keep == 1 || count <= log_size
    mat_neigh = mat_neigh * save_label;
    mat_neigh = double(mat_neigh > 0);
    count = count + 1;
    if mat_neigh == temp
        keep = 0;
    end
    temp = mat_neigh;
end
CV = conv_mat(mat_neigh,c);
fprintf('conv_mat function terminated \n');

newImg = tag_comp(mat_Tag, CV);
fprintf('tag_comp function terminated \n');

%disp(newImg);
imshow(newImg);
end

function [tag_mat, save_lab] = mat_transform(img)
%MAT_TRANSFORM this function run over all the pixels, and give a specificlabel to each potential connected component . It returns the "tagged matrix" and save_lab will be used to create the convertor vector later
 imgTemp = img;
[row,col] = size(img);
left = zeros(row+1,1);
up = zeros(1,col);
imgTemp = [up; imgTemp];
imgTemp = [left imgTemp];
labels = 1;
diag_label = ones(1,col*3);
save_label = diag(diag_label);
for r = 2:row + 1
        for c = 2:col + 1
            if imgTemp(r,c) ~= 0
                if imgTemp(r-1,c) ~= 0 
                    if imgTemp(r,c-1) ~= 0 
                        save_label(imgTemp(r-1,c),imgTemp(r,c-1)) = 1;
                        save_label(imgTemp(r,c-1), imgTemp(r-1,c)) = 1;
                    end
                    imgTemp(r,c) = imgTemp(r-1,c);
                elseif imgTemp(r,c-1) ~= 0
                    imgTemp(r,c) = imgTemp(r,c-1);
                else
                    imgTemp(r,c) = labels;
                    labels = labels + 1;
                end
            end
        end
end
lab = labels -1;
save_lab = save_label(1:lab, 1:lab);
%disp(save_lab);
tag_mat = imgTemp(2:row+1, 2:col+1);
%disp(tag_mat);
end

function [vec] = conv_mat(mat_neigh, n)
%CONV_MAT this function crete the converter vector from the labelled Marix
%in the way we learnt in the tirgul
CV = mat_neigh(1,:);
label = 2;

for i=2:n
    if CV(i) == 0
        CV = CV + mat_neigh(i,:) * label;
        %disp(CV);
        label = label + 1;
    end
end
vec = CV;
%disp(vec);
end

function [mat_lab] = tag_comp(mat_tag, CV)
%TAG_COMP this function take the converter vector and the  mat_tag matrix and
%create the sequenced connected component, i.e the final expected matrix.
imgTmp = mat_tag;
[row,col] = size(mat_tag);
for r = 1:row
        for c = 1:col
            if imgTmp(r,c) ~= 0
                imgTmp(r,c) = CV(imgTmp(r,c));
            end
        end
end
mat_lab = imgTmp;
%disp(mat_lab);
end

