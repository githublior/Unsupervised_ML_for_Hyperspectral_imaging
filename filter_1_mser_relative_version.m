function [f_1_r_regions] = filter_1_mser_relative_version(MserRegions,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist);
fprintf('first filter process: remove close aereas. \n');
indx_to_rm = [];
%s_loc = MserRegions.Count;
[s_loc, ~] = size(MserRegions.Location);
distance_matrix = zeros(s_loc,s_loc);
for i=1:s_loc
    fprintf('%d / %d \n',i, s_loc);
    if ismember(i, indx_to_rm) 
        continue
    end
    for j=i+1:s_loc
        if ismember(j, indx_to_rm)
            continue
        end
        
        loc_i= MserRegions(i).Location;
        nbr_px_i = MserRegions(i).Count;
        x_i =double(loc_i(1,1));
        y_i =double(loc_i(1,2));
        loc_j= MserRegions(j).Location;
        nbr_px_j = MserRegions(j).Count;
        x_j =double(loc_j(1,1));
        y_j =double(loc_j(1,2));
        distance_matrix(i,j)= abs(x_i - x_j) + abs(y_i - y_j);
        
        if distance_matrix(i,j)< min_dist && distance_matrix(i,j)>0
%        if (distance_matrix(i,j)<  (1-compress_percentage) * (sqrt(nbr_px_i/pi)+sqrt(nbr_px_j/pi)))  && distance_matrix(i,j)>0
            %if both same size (b/s): keep bigger
            % else: keep smaller
            
            ax_i = MserRegions(i).Axes;
            size_i =double(ax_i(1,1)) + double(ax_i(1,2));
            
            ax_j = MserRegions(j).Axes;
            size_j =double(ax_j(1,1)) + double(ax_j(1,2));
            % case same size
            ma = max(size_i,size_j);
            mi = min(size_i,size_j);
            
            if mi/ma  >= same_size_percentage
                %keep bigger
                if size_i >size_j
                    indx_to_rm =[indx_to_rm, j];
                else
                    indx_to_rm =[indx_to_rm, i];
                end
            % case one big one small    
            else
                %keep smaller
                if size_i >size_j
                    indx_to_rm =[indx_to_rm, i];
                else
                    indx_to_rm =[indx_to_rm, j];
                end
            end
            
            indx_to_rm =unique(indx_to_rm);
            
        end
        
    end
end
%{
indx_to_rm = [];
for i=1:s_loc
    for j=1:s_loc
        if i==j
            continue
        end
        if distance_matrix(i,j)<  (1-compress_percentage) * (sqrt(nbr_px_i/pi)+sqrt(nbr_px_j/pi))  && distance_matrix(i,j)>0
            %if both same size (b/s): keep bigger
            % else: keep smaller
            ax_i = MserRegions(i).Axes;
            size_i =double(ax_i(1,1)) + double(ax_i(1,2));
            
            ax_j = MserRegions(j).Axes;
            size_j =double(ax_j(1,1)) + double(ax_j(1,2));
            % case same size
            if size_i/size_j  >= same_size_percentage
                %keep bigger
                if size_i >size_j
                    indx_to_rm =[indx_to_rm, j];
                else
                    indx_to_rm =[indx_to_rm, i];
                end
            % case one big one small    
            else
                %keep smaller
                if size_i >size_j
                    indx_to_rm =[indx_to_rm, i];
                else
                    indx_to_rm =[indx_to_rm, j];
                end
            end
            
        end
    end
end
%}
size(indx_to_rm)
indx_to_rm =unique(indx_to_rm);
size(indx_to_rm)



%removing them
PixelList = [];
for i=1:size(MserRegions)
    if  ~ismember(i, indx_to_rm)
        PixelList = [PixelList ;{MserRegions(i).PixelList}];
    end
end

[~,nn]=size(indx_to_rm);
fprintf('%d objects removed with filter_1 \n', nn);

f_1_r_regions = MSERRegions(PixelList);


%PLOTTING
if plot_f_1
    %all_3_regr
    figure('Name','filter_1 ');imshow(bg); hold on;
    %figure; imshow(rgb); hold on;
    plot(f_1_r_regions);
%    plot(regions,'showPixelList',false,'showEllipses',true);

end