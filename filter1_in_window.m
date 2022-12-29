%more_index_to_rm
function  [PL,extracted] = filter1_in_window(PL,regions, local_indx_in_this_window, compress_percentage, same_size_percentage,min_dist);

%reduce mserRegions to only the one in local_indx_in_this_window
window_PixelList= [];
for indx= 1:size(local_indx_in_this_window)
    window_PixelList = [window_PixelList ; {regions(local_indx_in_this_window(indx)).PixelList}];
end
%class(window_PixelList)
%size(window_PixelList)
[this,~]=size(window_PixelList);
%fprintf('in this window we have %d regions to filter \n',this);
local_mserR = MSERRegions(window_PixelList);


%apply classic filter_1 to the local mserRegions.
more_index_to_rm = [];

[s_loc, ~] = size(local_indx_in_this_window);
distance_matrix = zeros(s_loc,s_loc);

for i=1:s_loc
    
    %print run time
    %fprintf('%d / %d \n',i, s_loc);
    
    if ismember(i, more_index_to_rm) 
        continue
    end
    for j=i+1:s_loc
        if ismember(j, more_index_to_rm)
            continue
        end
        
        loc_i= local_mserR(i).Location;
        [s_i,~] = size(local_mserR(i).PixelList);
        nbr_px_i = s_i;
        x_i =double(loc_i(1,1));
        y_i =double(loc_i(1,2));
        loc_j= local_mserR(j).Location;
        [s_j,~] = size(local_mserR(j).PixelList);
        nbr_px_j = s_j;
        x_j =double(loc_j(1,1));
        y_j =double(loc_j(1,2));
        distance_matrix(i,j)= abs(x_i - x_j) + abs(y_i - y_j);
        
        if distance_matrix(i,j)< min_dist && distance_matrix(i,j)>0
     %   if (distance_matrix(i,j)<  (1-compress_percentage) * (sqrt(nbr_px_i/pi)+sqrt(nbr_px_j/pi)))  && distance_matrix(i,j)>0
            %if both same size (b/s): keep bigger
            % else: keep smaller
            
            ax_i = local_mserR(i).Axes;
            size_i =double(ax_i(1,1)) + double(ax_i(1,2));
            
            ax_j = local_mserR(j).Axes;
            size_j =double(ax_j(1,1)) + double(ax_j(1,2));

            ma = max(size_i,size_j);
            mi = min(size_i,size_j);
            % case same size
            if mi/ma  >= same_size_percentage
                %keep bigger
                if size_i >size_j
                    more_index_to_rm =[more_index_to_rm, j];
%                    more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(j)];
                    PL(i)={unique([PL{(local_indx_in_this_window(i))} ; PL{(local_indx_in_this_window(j))} ] ,'rows')};


                else
                    more_index_to_rm =[more_index_to_rm, i];
 %                   more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(i)];
                    PL(j)={unique([PL{(local_indx_in_this_window(j))} ; PL{(local_indx_in_this_window(i))} ] ,'rows')};


                end
            % case one big one small    
            else
                %keep smaller
                if size_i >size_j
                    more_index_to_rm =[more_index_to_rm, i];
%                    more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(i)];
%                    PL(j)={unique([PL{(local_indx_in_this_window(j))} ; PL{(local_indx_in_this_window(i))} ] ,'rows')};

                    
                else
                    more_index_to_rm =[more_index_to_rm, j];
%                    more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(j)];
%                    PL(i)={unique([PL{(local_indx_in_this_window(i))} ; PL{(local_indx_in_this_window(j))} ] ,'rows')};

                end
            end
            
            more_index_to_rm =unique(more_index_to_rm);

        end
        
    end
end



[~,w]= size(more_index_to_rm);

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


extracted = [];
for indx=1:w
    extracted = [extracted;local_indx_in_this_window(more_index_to_rm(indx)) ];
end
%}
fprintf(' %d new elt to remove in this window \n', w);


%{
%  ? CAN REMOVE ? 
fprintf('rm \n');
size(more_index_to_rm)
more_index_to_rm =unique(more_index_to_rm);
size(more_index_to_rm)
fprintf('rm \n');
%}