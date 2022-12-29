
function [f_1_r_regions] = filter_1_mser_v3(MserRegions,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,a,b,max_area_raw, max_area_col,min_range,max_range );
%{
this function filter the mser region. removes all regions that are included
in one another .
in this version, in order to have an acceptable run time, we want to work over a sliding
window on the wkn_img of (2*max_aerea) * (2*max_aerea) moving by a stride
of max_aerea. (hozitontal and vertical.)
%}

indx_to_rm = [];
PL = MserRegions.PixelList;

%order regions locations x and  y while keeping their original row index
indexed_loc=[];
for temp1=1:MserRegions.Count
    loc_xandy = MserRegions(temp1).Location;
    x_l =double(loc_xandy(1,1));
    y_l =double(loc_xandy(1,2));
    indexed_loc = [indexed_loc; [temp1 x_l y_l]];
end

x_ordered_loc  = sortrows(indexed_loc,2);
y_ordered_loc  = sortrows(indexed_loc,3);



start_raw = 1;
start_col = 1;

%{
for raw_wind=max_area:max_aera/2:a
    for col_wind=max_aera:max_aera/2:b
        local_indx_to_treat = [];

        if x_ordered_loc(raw_check,2)> raw_wind
            break
        else if x_ordered_loc(raw_check,2)<start_raw
                continue
             else
                local_indx_to_treat = [local_indx_to_treat; MserRegions().Location  ]
            end
        end
        for raw_check=x_indx:MserRegions.Count

        
        end
    
        
        start_col= col_wind;
    end
    
    start_raw= raw_wind;
end
 %}

%slide window by row
for raw_wind=min(a,max_area_raw):max_area_raw/2:a
    if a-raw_wind < max_area_raw/2
        raw_wind = a;
    end
    y_indx_window=[];
    % check if x into this row_window for all elmt of mserregions
    for raw_check=1:MserRegions.Count
        if y_ordered_loc(raw_check,3)> raw_wind
            break
        else
            if y_ordered_loc(raw_check,3)<start_raw
                continue
            else
                
                y_indx_window= [y_indx_window; y_ordered_loc(raw_check,1)];
            end
        end
    end
        

    %slide window by col
    for col_wind=min(b,max_area_col):max_area_col/2:b
        if b-col_wind < max_area_col/2
            col_wind = b;
        end
        x_indx_window=[];
        % check if y into this row_window for all elmt of mserregions
        for col_check=1:MserRegions.Count
            if x_ordered_loc(raw_check,2)> col_wind
                break
            else
                if x_ordered_loc(raw_check,2)< start_col
                    continue
                else
                    x_indx_window = [x_indx_window;  x_ordered_loc(col_check,1)];
                    
                end
            end
        end
        

        %select index loc that match x and y _indx_wondow
        [local_indx_in_this_window,~]=intersect(x_indx_window,y_indx_window);
        
        [this,~] = size(local_indx_in_this_window);
        fprintf(' %d   potential regions in the window %d , %d, %d , %d \n', this, start_raw , raw_wind ,start_col ,  col_wind);
        
        if this >0         

            % apply the filter 
            [PL, more_index_to_rm] = filter1_in_window(PL,MserRegions,local_indx_in_this_window, compress_percentage, same_size_percentage,min_dist);
%            more_index_to_rm = filter1_in_window(PL,MserRegions,local_indx_in_this_window, compress_percentage, same_size_percentage,min_dist);
            [first , ~] =size(more_index_to_rm);

            for elt=1:first
                indx_to_rm = [indx_to_rm,more_index_to_rm(elt) ];
            end        
            fprintf('before');
            size(indx_to_rm)
            indx_to_rm = unique(indx_to_rm);
            fprintf('after');
            size(indx_to_rm)
        end
        start_col= col_wind;
    end
    
    
    start_raw= raw_wind;
end
        



[t,y] = size(PL);

fprintf('%d %d PL size' , t,y);

%removing them
PixelList = [];
for i=1:size(MserRegions)
    if  ~ismember(i, indx_to_rm)
        PixelList = [PixelList ;PL(i)];
%        PixelList = [PixelList ;{MserRegions(i).PixelList}];
    end
end

class(PixelList)
size(PixelList)

[~,nn]=size(indx_to_rm);
fprintf('%d objects removed with filter_1 \n', nn);

f_1_r_regions = MSERRegions(PixelList);

%f_1_r_regions = detectMSERFeatures(PixelList, 'RegionAreaRange',[min_range ,max_range]);




%PLOTTING
if plot_f_1
    %all_3_regr
    figure();imshow(bg); hold on;
    hp = impixelinfo;
    %figure; imshow(rgb); hold on;
    %plot(f_1_r_regions);
    %    plot(f_1_r_regions,'showPixelList',true);
    plot(f_1_r_regions,'showPixelList',true,'showEllipses',true);

end