% image process from kmeans 10 output
%given an input image taken from kmeans algorithm, we want to detect cels
% using connected component. 
% we choose the cluster representing the cells kernel, 
% apply tag connected component function
% filter the output with kernel max size - min size
% iteratively check the presence of the cytoplasm's cluster 



%get output of previous run of this code 
%kernel_out = imread('kernel_out.png');


%filtered_selection
%{
%parameters
% SELECT working_img : which kmeans to work on . 
%{
%temp = imread('10_kmeans_of_mat_P055A_X2000Y0_R1.png');
%crop_k =mat2gray(temp);

%first img loads
temp = imread('10K_means_out.png');
cropk= imcrop(temp,[520 400  (819-520) (699-400)]);
%imtool(cropk);
%}

%reference cell
wkn_img = imread('C161VA_x20_X700Y600_20_kmeans.png');
%crop_k=imcrop(wkn_img,[520 400  (819-520) (699-400)]);
name = 'C161VA_x20_X700Y600_20_kmeans';

% select  clusters to filter 
%{
%first img clusters values
first_kernel = [0,28,142,198,170];
first_cytoplasm = [227,255];
first_membrane = [55,113,57];
%}

%reference cell clusters values

%C161VA_x20_X700Y600_20_kmeans
k_13 =[13];
reference_kernel = [13,121, 242, 40];
reference_cyto = [201,161, 148];
c_201_242 = [201,242];
c_40 = [40];
c_134 = [134];
c_40_242 = [40,242];
c_148 = [148];
temp  =[ 40, 242];
cancerous_zone = [40,134 ,201 , 148,242];
canc_cell =[.16];
cell_area = c_40_242;
cell_area_name = 'c_40_242';


% select min - max allowed size
Max_allowed_cc_size = 150000000;
Min_allowed_cc_size = 200;


%select background
%{
%kernel_out_8_300 = imread('kernel_out_8_300.png');
%support_name= 'kernel_8_300';
%support=mat2gray(kernel_out_8_300);

%m_k_8_300 = imread('kernel_8_300_on_empty_filtered_img.png');
%support_name= 'm_k_8_300';
%support=mat2gray(m_k_8_300);

%membrane_kernel_8_100 = imread('kernel_out_8_100.png');
%support =mat2gray(membrane_kernel_8_100);
%support_name= 'membrane_kernel_8_100';


[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;


background = mat2gray(imread('C161VA_x20_X700Y600_20_kmeans_k_13_20_1000_on_empty.png'));
bg_name= 'k13_20_1000';

%}
[a,b] = size(wkn_img);
binay =zeros(a,b);
bg_name= 'empty';
background = binay;



bg_reduction = 0;
save_filtered_img = 0;


final_name = append(name, '_', cell_area_name,'_',int2str(Min_allowed_cc_size),'_',int2str(Max_allowed_cc_size),'_on_',bg_name);

f_i = filtered_selection( wkn_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area,save_filtered_img,final_name,bg_reduction);
imtool(f_i);

%{
kernel = -1*ones(3);
kernel(2,2) = 17;
enhancedImage = imfilter(f_i, kernel);
enhancedImage = imfill(enhancedImage,'holes');
imtool(enhancedImage);
imwrite(enhancedImage,'C161VA_x20_X700Y600_cancerous_zone.png');
%}

%{

function [filtered_img] = temp_filtered_selection(working_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area, save_filtered_img, name,bg_reduction)
%{
working_img : kmeans mat we work on
background : if the background is empty or contains previous filters.
Min_allow, Max_allow_cc : treshold size of connected components
cell_area: slected clusters list
 save_filtered_img: 1 to save img , 0 else
 name: is save img, name of it
bg_reduction: if the BG contain previous filters, 1 if we want to give them another color, 0 else
%}
 Min_allowed_cc_size
 Max_allowed_cc_size
 cell_area
 save_filtered_img
 name,bg_reduction
    km = imread('C161VA_x20_X700Y600_20_kmeans.png');
    [tr1,tr2] =size(background);
    
    %reduce supports value
    if bg_reduction ==1
        for qwe=1:tr1
            for wer=1:tr2
                 background(qwe,wer) = max( background(qwe,wer) - .3,0);
            end
        end
    end
    %imtool(support);
    

    [a,b] = size(working_img);
    if tr1 ~= a || tr2 ~=b
        fprintf('dim err: support matrix have wrong dimensions \n');
    end

    %convert image to binary according to chosens cluster

    %filter the relevant cluster : set their value to 1, other pixel 0
    binay=zeros(a,b);
    binay = select_specific_clusters( binay, working_img, cell_area);


    %PRINTING AND WRITING ZONE
%{
    mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/3.rawstack.mat');
    %imshow(mymat.RGB);
    crop_RGB = imcrop(mymat.RGB,[520 400  (819-520) (699-400)]);
    rgb= imtool(crop_RGB);
    set(rgb,'NumberTitle','off','Name','rgb');

    cropk  = imtool(crop_k);
    set(cropk,'NumberTitle','off','Name','crop_k');

    %ko = imtool(kernel_out);
    %set(ko,'NumberTitle','off','Name','kernel out');
%}
    bin = imtool(mat2gray(binay));
    set(bin,'NumberTitle','off','Name','any size clusters mask ');

    
    %'all cyto+ filtered kernel

    %imtool(mat2gray(kernel_out));
    %imwrite(binay, 'binay1.png');

    

    %APPLY CONNECTED COMPONENT FUNCTIONS.
    %https://fr.mathworks.com/help/images/ref/regionprops.html#d123e240356
    %https://fr.mathworks.com/help/images/ref/bwconncomp.html#d123e31653

    [L,~] = bwlabel(binay);
    %imtool(mat2gray(L));

    % count # pixels in each CC of L.
    map = nbr_elt_in_clusters(a,b,L);
    %{
    keymat = cell2mat(keys(map));
    freq = cell2mat(values(map));
    [~, sm]=size(keymat);
    %}
    cinay=zeros(a,b);
    cinay = select_specific_clusters( cinay, working_img, [40]);



    CC = bwconncomp(cinay);
    %CC.NumObjects
    %CC.PixelIdxList
    stats = regionprops('Table',CC,'basic');
    centroid_list = stats.Centroid;
    [one, ~] = size(centroid_list);
   %{ 
ctrd = centroid_list(1,:);

    BoundingBox = stats.BoundingBox;
    fprintf('bb %f \n', BoundingBox(1,:));
    [ca, ~] = size(centroid_list);

    %}
    
    % select connected component with appropriate size, and output them in "filtered_img"
    %filtered_img=zeros(a,b);
    filtered_img=background;
    for un= 1:a
        for deux =1:b
            if isKey(map, L(un,deux))
                if map(L(un,deux)) <Max_allowed_cc_size  && map(L(un,deux))> Min_allowed_cc_size
                    bbol = 0;
                    for mtemp =1: one
                    if abs(un- centroid_list(mtemp,1)) +abs(deux- centroid_list(mtemp,2))  <50
                        bbol=1;
                    end
                    end
                    if bbol ==0
                        filtered_img(un,deux)=1;
                    end
                end
            end
            %{
                %if cell2mat(values(map,L(un,deux))) <Max_allowed_cc_size  && freq( L(un,deux))>Min_allowed_cc_size
             %   filtered_img(un,deux)=1;
            %end
            %fprintf( L(un,deux) )
            %for trois=1:sm
            %    if L(un,deux) ==  keymat(trois)
            %        if freq(trois) <Max_allowed_cc_size  && freq(trois)>Min_allowed_cc_size
             %           filtered_img(un,deux)=1;
             %       end
             %end
             %}
         end
     end

    %{

    CC = bwconncomp(binay);
    %CC
    %CC.PixelIdxList
    stats = regionprops('Table',CC,'basic');
    centroid_list = stats.Centroid;
    number_of_elt_in_ccs = stats.Area;
    BoundingBox = stats.BoundingBox;
    [ca, ~] = size(centroid_list);

    count= 0;
    for i = 1:CC.NumObjects
        if number_of_elt_in_ccs(i,1) < Max_allowed_cc_size && number_of_elt_in_ccs(i,1) > Min_allowed_cc_size

            count = count +1;
            %fprintf('connected compo idx : %d \n',i);
            ctrd = centroid_list(i,:);
            ctrd_x = round(ctrd(1));
            ctrd_y = round(ctrd(2));
            conec_compo_index= L(ctrd_x,ctrd_y);
        %if conec_compo_index == 155
            %{
            fprintf('nbr elt %f \n',number_of_elt_in_ccs(i));
            fprintf('centroid: %f \n', centroid_list(i,:));
            %fprintf('bb %f \n', BoundingBox(i,:));
            %}
            x = round(BoundingBox(i,1));
            y= round(BoundingBox(i,2));
            W = BoundingBox(i,3);
            H = BoundingBox(i,4);
           % if ctr1_x < 224 && ctr1_x > 214 && ctr1_y > 160 && ctr1_y < 170
                        %craft
            %fprintf('x %d, y %d, W %d, H %d \n', x,y,W,H);

            L_vals_in_cc = get_L_vals_and_occurences_mapping(a,b,y,x,W,H,L);  
            %BIZARRE x,y inverser

            keymat = cell2mat(keys(L_vals_in_cc));
            freq = cell2mat(values(L_vals_in_cc));
            [~,l2] = size(freq);
            %remove(L_vals_in_cc,0);


           % detect wich key has the highest occurence by checking the
           % index o the highest value
           [maxfreq, maxindx] = get_maxfreq_max_indx(keymat,freq,l2);
           if l2 >0
               %fprintf('maxindx %d , maxval %d \n', maxindx, maxfreq);
               %fprintf('chosen L value: %d  \n', keymat(maxindx));
               correct_L_label =  keymat(maxindx);
               if correct_L_label == 154
                   fprintf('x %d, y %d, W %d, H %d \n', x,y,W,H);
                   fprintf('%d \n' ,number_of_elt_in_ccs(i,:));
               end
               for j=1:a
                   %j = x : x+W
                for  k=1:b
                    %k = y : y+H
                    t1=min(j,a);
                    t2 = min(k,b);
                    if L(t2,t1) == correct_L_label
                        %BIZARRE t2,t1 inverser
                        filtered_img(t2,t1)=1;
                         %BIZARRE t2,t1 inverser
                    end
                end
               end
           end
        end
            %fprintf('index %d, centro %d, %d , concompo_ind_L: %d and value on L %f \n',i,ctr1_x , ctr1_y  , conec_compo_index, L(181,294));

    %{
            for j = x : x+W
                for k = y : y+H
                    t1=min(j,a);
                    t2 = min(k,b);

                    if binay(t1,t2) == 0 && L(t1,t2)~= 0
                        fprintf('0 bin vs  ~0 L missmatch : %d %d \n', t1,t2);
                        countmiss = countmiss +1;
                    else
                        if binay(t1,t2) ~= 0 && L(t1,t2)== 0
                            fprintf('1 bin vs 0 L missmatch : %d %d \n', t1,t2);
                            countmiss = countmiss +1;
                         end   
                    end

                    %filtered_img(t1,t2) = 1;
                    filtered_img(t1,t2) = binay(t1,t2);
                end
            end 
    %}
       % end
    end
    %}
    
    %{
    BWdfill = imfill(filtered_img,'holes');
    imtool(BWdfill);
   
    %add treshold: remove every cc that doestn fit minimum cell size.
    [L2,~] = bwlabel(BWdfill);
    imtool(L2);
    % count # pixels in each CC of L.
    second_map = nbr_elt_in_clusters(a,b,L2);
   
    keymat2 = cell2mat(keys(second_map));
    freq2 = cell2mat(values(second_map));
    [~, sm2]=size(keymat2);
    Min_cell_size= 150;
    new_filtered_img= filtered_img;
    for un= 1:a
        for deux =1:b
            for trois=1:sm2
                
                if L(un,deux) ==  keymat2(trois)
                    fprintf('%d, x= %d , y= %d \n',freq2(trois),un,deux );
                    if freq2(trois) <Min_cell_size 
                        new_filtered_img(un,deux)=0;
                    end
                end
            end
        end
    end
    

    imtool(new_filtered_img)
    %}
    if save_filtered_img==1
        output_name =append(name,'.png' );
        imwrite(filtered_img, output_name);
    end
    
    fprintf('filter selection function terminated \n');
end


%fprintf('nuber of accpeted connected component : %d \n',count);
%hFigure = imtool(mat2gray(filtered_img));
%set(hFigure,'NumberTitle','off','Name','filtered cytoplasm U kernel');
%imwrite(filtered_img, 'kernel_out.png')


function [bg] = select_specific_clusters( bg, working_img, cell_area)
    %{
    %for true cell analysis
    blank: .22, .56, .67 little .33
    cyto : .78 1 , .89
    kernel/membrane : .44
    between kernel-cyto :0 
    %}

    % miss between kernel cytoplasm : 57 - 85
    [a,b] = size(working_img);
    [~,c] = size(cell_area);

    for i = 1:a
        for j= 1:b
            bol =0;
            for k= 1:c
                if uint8(working_img(i,j)) == cell_area(k)
                    bol = 1;
                end
            end
            if bol ==1
                bg(i,j)=1;
            end
            %else
             %   binay(i,j) = 0; 
                %binay(i,j) = mat2gray(kernel_out(i,j));
                %binay(i,j) = max(mat2gray(kernel_out(i,j))-.5,0);
        end
    end
end

function [ L_map ] = nbr_elt_in_clusters(a,b, L) 
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
end

%}

%}
function [filtered_img] = filtered_selection(working_img,background, Min_allowed_cc_size, Max_allowed_cc_size,cell_area, save_filtered_img, name,bg_reduction)
%{
working_img : kmeans mat we work on
background : if the background is empty or contains previous filters.
Min_allow, Max_allow_cc : treshold size of connected components
cell_area: slected clusters list
 save_filtered_img: 1 to save img , 0 else
 name: is save img, name of it
bg_reduction: if the BG contain previous filters, 1 if we want to give them another color, 0 else
%}
    [tr1,tr2] =size(background);
    
    %reduce supports value
    if bg_reduction ==1
        for qwe=1:tr1
            for wer=1:tr2
                 background(qwe,wer) = max( background(qwe,wer) - .3,0);
            end
        end
    end
    %imtool(support);
    

    [a,b] = size(working_img);
    if tr1 ~= a || tr2 ~=b
        fprintf('dim err: support matrix have wrong dimensions \n');
    end

    %convert image to binary according to chosens cluster

    %filter the relevant cluster : set their value to 1, other pixel 0
    binay=zeros(a,b);
    binay = select_specific_clusters( binay, working_img, cell_area);


    %PRINTING AND WRITING ZONE
%{
    mymat = load('/Users/lior/Desktop/wetransfer_3-rawstack-jpg_2021-08-30_1417/3.rawstack.mat');
    %imshow(mymat.RGB);
    crop_RGB = imcrop(mymat.RGB,[520 400  (819-520) (699-400)]);
    rgb= imtool(crop_RGB);
    set(rgb,'NumberTitle','off','Name','rgb');

    cropk  = imtool(crop_k);
    set(cropk,'NumberTitle','off','Name','crop_k');

    %ko = imtool(kernel_out);
    %set(ko,'NumberTitle','off','Name','kernel out');
%}
    bin = imtool(mat2gray(binay));
    set(bin,'NumberTitle','off','Name','any size clusters mask ');

    
    %'all cyto+ filtered kernel

    %imtool(mat2gray(kernel_out));
    %imwrite(binay, 'binay1.png');

    

    %APPLY CONNECTED COMPONENT FUNCTIONS.
    %https://fr.mathworks.com/help/images/ref/regionprops.html#d123e240356
    %https://fr.mathworks.com/help/images/ref/bwconncomp.html#d123e31653

    [L,~] = bwlabel(binay);
    %imtool(mat2gray(L));

    % count # pixels in each CC of L.
    map = nbr_elt_in_clusters(a,b,L);
    %{
    keymat = cell2mat(keys(map));
    freq = cell2mat(values(map));
    [~, sm]=size(keymat);
    %}
    
    % select connected component with appropriate size, and output them in "filtered_img"
    %filtered_img=zeros(a,b);
    filtered_img=background;
    for un= 1:a
        for deux =1:b
            if isKey(map, L(un,deux))
                if map(L(un,deux)) <Max_allowed_cc_size  && map(L(un,deux))> Min_allowed_cc_size
                    filtered_img(un,deux)=1;
                end
            end
            %{
                %if cell2mat(values(map,L(un,deux))) <Max_allowed_cc_size  && freq( L(un,deux))>Min_allowed_cc_size
             %   filtered_img(un,deux)=1;
            %end
            %fprintf( L(un,deux) )
            %for trois=1:sm
            %    if L(un,deux) ==  keymat(trois)
            %        if freq(trois) <Max_allowed_cc_size  && freq(trois)>Min_allowed_cc_size
             %           filtered_img(un,deux)=1;
             %       end
             %end
             %}
         end
     end

    %{

    CC = bwconncomp(binay);
    %CC
    %CC.PixelIdxList
    stats = regionprops('Table',CC,'basic');
    centroid_list = stats.Centroid;
    number_of_elt_in_ccs = stats.Area;
    BoundingBox = stats.BoundingBox;
    [ca, ~] = size(centroid_list);

    count= 0;
    for i = 1:CC.NumObjects
        if number_of_elt_in_ccs(i,1) < Max_allowed_cc_size && number_of_elt_in_ccs(i,1) > Min_allowed_cc_size

            count = count +1;
            %fprintf('connected compo idx : %d \n',i);
            ctrd = centroid_list(i,:);
            ctrd_x = round(ctrd(1));
            ctrd_y = round(ctrd(2));
            conec_compo_index= L(ctrd_x,ctrd_y);
        %if conec_compo_index == 155
            %{
            fprintf('nbr elt %f \n',number_of_elt_in_ccs(i));
            fprintf('centroid: %f \n', centroid_list(i,:));
            %fprintf('bb %f \n', BoundingBox(i,:));
            %}
            x = round(BoundingBox(i,1));
            y= round(BoundingBox(i,2));
            W = BoundingBox(i,3);
            H = BoundingBox(i,4);
           % if ctr1_x < 224 && ctr1_x > 214 && ctr1_y > 160 && ctr1_y < 170
                        %craft
            %fprintf('x %d, y %d, W %d, H %d \n', x,y,W,H);

            L_vals_in_cc = get_L_vals_and_occurences_mapping(a,b,y,x,W,H,L);  
            %BIZARRE x,y inverser

            keymat = cell2mat(keys(L_vals_in_cc));
            freq = cell2mat(values(L_vals_in_cc));
            [~,l2] = size(freq);
            %remove(L_vals_in_cc,0);


           % detect wich key has the highest occurence by checking the
           % index o the highest value
           [maxfreq, maxindx] = get_maxfreq_max_indx(keymat,freq,l2);
           if l2 >0
               %fprintf('maxindx %d , maxval %d \n', maxindx, maxfreq);
               %fprintf('chosen L value: %d  \n', keymat(maxindx));
               correct_L_label =  keymat(maxindx);
               if correct_L_label == 154
                   fprintf('x %d, y %d, W %d, H %d \n', x,y,W,H);
                   fprintf('%d \n' ,number_of_elt_in_ccs(i,:));
               end
               for j=1:a
                   %j = x : x+W
                for  k=1:b
                    %k = y : y+H
                    t1=min(j,a);
                    t2 = min(k,b);
                    if L(t2,t1) == correct_L_label
                        %BIZARRE t2,t1 inverser
                        filtered_img(t2,t1)=1;
                         %BIZARRE t2,t1 inverser
                    end
                end
               end
           end
        end
            %fprintf('index %d, centro %d, %d , concompo_ind_L: %d and value on L %f \n',i,ctr1_x , ctr1_y  , conec_compo_index, L(181,294));

    %{
            for j = x : x+W
                for k = y : y+H
                    t1=min(j,a);
                    t2 = min(k,b);

                    if binay(t1,t2) == 0 && L(t1,t2)~= 0
                        fprintf('0 bin vs  ~0 L missmatch : %d %d \n', t1,t2);
                        countmiss = countmiss +1;
                    else
                        if binay(t1,t2) ~= 0 && L(t1,t2)== 0
                            fprintf('1 bin vs 0 L missmatch : %d %d \n', t1,t2);
                            countmiss = countmiss +1;
                         end   
                    end

                    %filtered_img(t1,t2) = 1;
                    filtered_img(t1,t2) = binay(t1,t2);
                end
            end 
    %}
       % end
    end
    %}
    
    %{
    BWdfill = imfill(filtered_img,'holes');
    imtool(BWdfill);
   
    %add treshold: remove every cc that doestn fit minimum cell size.
    [L2,~] = bwlabel(BWdfill);
    imtool(L2);
    % count # pixels in each CC of L.
    second_map = nbr_elt_in_clusters(a,b,L2);
   
    keymat2 = cell2mat(keys(second_map));
    freq2 = cell2mat(values(second_map));
    [~, sm2]=size(keymat2);
    Min_cell_size= 150;
    new_filtered_img= filtered_img;
    for un= 1:a
        for deux =1:b
            for trois=1:sm2
                
                if L(un,deux) ==  keymat2(trois)
                    fprintf('%d, x= %d , y= %d \n',freq2(trois),un,deux );
                    if freq2(trois) <Min_cell_size 
                        new_filtered_img(un,deux)=0;
                    end
                end
            end
        end
    end
    

    imtool(new_filtered_img)
    %}
    if save_filtered_img==1
        output_name =append(name,'.png' );
        imwrite(filtered_img, output_name);
    end
    
    fprintf('filter selection function terminated \n');
end


%fprintf('nuber of accpeted connected component : %d \n',count);
%hFigure = imtool(mat2gray(filtered_img));
%set(hFigure,'NumberTitle','off','Name','filtered cytoplasm U kernel');
%imwrite(filtered_img, 'kernel_out.png')


function [bg] = select_specific_clusters( bg, working_img, cell_area)
    %{
    %for true cell analysis
    blank: .22, .56, .67 little .33
    cyto : .78 1 , .89
    kernel/membrane : .44
    between kernel-cyto :0 
    %}

    % miss between kernel cytoplasm : 57 - 85
    [a,b] = size(working_img);
    [~,c] = size(cell_area);

    for i = 1:a
        for j= 1:b
            bol =0;
            for k= 1:c
                if uint8(working_img(i,j)) == cell_area(k)
                    bol = 1;
                end
            end
            if bol ==1
                bg(i,j)=1;
            end
            %else
             %   binay(i,j) = 0; 
                %binay(i,j) = mat2gray(kernel_out(i,j));
                %binay(i,j) = max(mat2gray(kernel_out(i,j))-.5,0);
        end
    end
end

function [ L_map ] = nbr_elt_in_clusters(a,b, L) 
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
end