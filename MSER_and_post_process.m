
function [f_2_regions] = MSER_and_post_process(wkn_img,min_range , max_range, f_1_min_distance,f_2_treshold_miss,plot_mser, plot_f_1, plot_f_2, bg);

    %MSER
    %[regions,cc] cc maybe usefull
    [regions,~] = detectMSERFeatures(wkn_img, 'RegionAreaRange',[min_range ,max_range]);
    fprintf('%d object detect at the begining of MSER \n', regions.Count);

    %PLOTTING mser
    if plot_mser
        figure('Name','MSER on chosen bg');imshow(bg); hold on;
        %figure('Name','MSER on RGB','NumberTitle','off');imshow(rgb); hold on;
        %plot(regions);
        %plot(regions(1),'showPixelList', true);
        plot(regions,'showEllipses',true);
        %plot(regions,'showPixelList',true,'showEllipses',true);
    end
    



    %POST PROCESS

    %filter 1: when some centroide are close , keep the biggest ellipse only

    f_1_min_distance =25;

    [s_loc, ~] = size(regions.Location);
    distance_matrix = zeros(s_loc,s_loc);
    for i=1:s_loc
        for j=1:s_loc
            loc_i= regions(i).Location;
            x_i =double(loc_i(1,1));
            y_i =double(loc_i(1,2));
            loc_j= regions(j).Location;
            x_j =double(loc_j(1,1));
            y_j =double(loc_j(1,2));
            distance_matrix(i,j)= abs(x_i - x_j) + abs(y_i - y_j);

        end
    end


    indx_to_rm_f_1 = [];
    for i=1:s_loc
        for j=1:s_loc
            if i==j
                continue
            end
            if distance_matrix(i,j)< f_1_min_distance
                ax_i = regions(i).Axes;
                size_i =double(ax_i(1,1)) + double(ax_i(1,2));

                ax_j = regions(j).Axes;
                size_j =double(ax_j(1,1)) + double(ax_j(1,2));

                if size_i >size_j
                    indx_to_rm_f_1 =[indx_to_rm_f_1, j];
                else
                    indx_to_rm_f_1 =[indx_to_rm_f_1, i];

                end

            end
        end
    end
    indx_to_rm_f_1 =unique(indx_to_rm_f_1);
    size(indx_to_rm_f_1);

    %removing them
    PixelList = [];
    for i=1:size(regions)
        if  ~ismember(i, indx_to_rm_f_1)
            PixelList = [PixelList ;{regions(i).PixelList}];
        end
    end

    [~,nn]=size(indx_to_rm_f_1);

    f_1_regions = MSERRegions(PixelList);
    fprintf('%d object detect  after filter_1  (%d objects removed) \n', f_1_regions.Count, nn);

    %PLOTTING filter_1
    if plot_f_1
        %all_3_regr
        figure('Name','filter_1 on chosen bg');imshow(bg); hold on;
        %figure; imshow(rgb); hold on;
        plot(f_1_regions);
        %plot(re_r,'showPixelList', true);
    end



    %filter_2: removing ellipses that surround void

    % percentage of miss to reject the element
    treshold_miss=f_2_treshold_miss ;

    indx_to_rm = [];
    nbr_elt = f_1_regions.Count;
    for i=1:nbr_elt
        [pxl_grp,~] = size(f_1_regions.PixelList(i));
        %fprintf('new group %d  \t there are %d pixels in this group \n', i,pxl_grp);

        miss_pixel=0;
        for j=1:pxl_grp
            %fprintf('new pixel %d \n',j);
            pxl_of_elt = f_1_regions.PixelList(i);
            y=(pxl_of_elt(j,1));
            x=(pxl_of_elt(j,2));
            if wkn_img(x,y) == 0
                miss_pixel = miss_pixel + 1 ;
            end
        end
        if double(miss_pixel)/pxl_grp > treshold_miss
            indx_to_rm =[indx_to_rm, i];
        end
        %size(indx_to_rm)
    end

    %removing them
    PixelList = [];
    for i=1:nbr_elt
        if  ~ismember(i, indx_to_rm)
            PixelList = [PixelList ;{f_1_regions(i).PixelList}];
        end
    end
    [~,nn]=size(indx_to_rm);

    f_2_regions = MSERRegions(PixelList);
    fprintf('%d object detect  after filter_2  (%d objects removed) \n', f_2_regions.Count, nn);



    %PLOTTING f_2
    if plot_f_2
        %all_3_regr
        figure('Name','filter_2 on chosen bg');imshow(bg); hold on;
        %figure; imshow(all_3_regr); hold on;
        plot(f_2_regions);
        plot(f_2_regions,'showPixelList', true);
    end

    fprintf('function MSER_and_post_process terminated');

