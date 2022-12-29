function [f_1_regions] = filter_1_mser(MserRegions, min_dist,plot_f_1,bg);

minimum_distance =min_dist;

[s_loc, ~] = size(MserRegions.Location);
distance_matrix = zeros(s_loc,s_loc);
for i=1:s_loc
    for j=1:s_loc
        loc_i= MserRegions(i).Location;
        x_i =double(loc_i(1,1));
        y_i =double(loc_i(1,2));
        loc_j= MserRegions(j).Location;
        x_j =double(loc_j(1,1));
        y_j =double(loc_j(1,2));
        distance_matrix(i,j)= abs(x_i - x_j) + abs(y_i - y_j);
        
    end
end


indx_to_rm = [];
for i=1:s_loc
    for j=1:s_loc
        if i==j
            continue
        end
        if distance_matrix(i,j)< minimum_distance
            ax_i = MserRegions(i).Axes;
            size_i =double(ax_i(1,1)) + double(ax_i(1,2));
            
            ax_j = MserRegions(j).Axes;
            size_j =double(ax_j(1,1)) + double(ax_j(1,2));
            
            if size_i >size_j
                indx_to_rm =[indx_to_rm, j];
            else
                indx_to_rm =[indx_to_rm, i];

            end
            
        end
    end
end
indx_to_rm =unique(indx_to_rm);
size(indx_to_rm);

%removing them
PixelList = [];
for i=1:size(MserRegions)
    if  ~ismember(i, indx_to_rm)
        PixelList = [PixelList ;{MserRegions(i).PixelList}];
    end
end

[~,nn]=size(indx_to_rm);
fprintf('%d objects removed with filter_1 \n', nn);

f_1_regions = MSERRegions(PixelList);


%PLOTTING
if plot_f_1
    %all_3_regr
    figure('Name','filter_1 ');imshow(bg); hold on;
    %figure; imshow(rgb); hold on;
    plot(f_1_regions);
    %plot(re_r,'showPixelList', true);
end