
function [W_PL , f_1_r_regions] = filter_doublon_intra(MserRegions,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,area_size_row,area_size_col,slide_by_raw, slide_by_col,min_range,max_range,close_enough);
%{
this function filter the mser region. removes all regions that are included
in one another .
in this version, in order to have an acceptable run time, we want to work over a sliding
window on the wkn_img of (2*max_aerea) * (2*max_aerea) moving by a stride
of max_aerea. (hozitontal and vertical.)
%}
conteur = 0 ;
indx_to_rm = [];
W_PL = [];

%order regions locations x and  y while keeping their original row index
indexed_loc=[];
for region_index=1:MserRegions.Count
    loc_xandy = MserRegions(region_index).Location;
    x_l =double(loc_xandy(1,1));
    y_l =double(loc_xandy(1,2));
    indexed_loc = [indexed_loc; [region_index x_l y_l]];
end

x_ordered_loc  = sortrows(indexed_loc,2);
y_ordered_loc  = sortrows(indexed_loc,3);




%slide window by row
start_raw = 1;
for raw_wind=min(area_size_row, slide_by_raw): slide_by_raw / 2 : area_size_row
    if area_size_row - raw_wind < slide_by_raw / 2
        raw_wind = area_size_row;
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
        
    start_col = 1;
    %slide window by col
    for col_wind=min(area_size_col,slide_by_col):slide_by_col/2:area_size_col
        if area_size_col-col_wind < slide_by_col/2
            col_wind = area_size_col;
        end
        x_indx_window=[];
        % check if y into this row_window for all elmt of mserregions
        for col_check=1:MserRegions.Count
            if x_ordered_loc(col_check,2)> col_wind
                break
            else
                if x_ordered_loc(col_check,2)< start_col
                    continue
                else
                    x_indx_window = [x_indx_window;  x_ordered_loc(col_check,1)];
                    
                end
            end
        end
        
        
        %vizualize where do we currently work:
        %imshow(bg);
        %r = drawrectangle('Position', [start_raw start_col (raw_wind-start_raw) (col_wind-start_col)],'StripeColor','r');
 

        %select index loc that match x and y _indx_wondow
        [local_indx_in_this_window,~]=intersect(x_indx_window,y_indx_window);

        [nbr_regions_in_window,~] = size(local_indx_in_this_window);
       % fprintf(' %d   potential regions in the window %d , %d, %d , %d \n', nbr_regions_in_window, start_raw , raw_wind ,start_col ,  col_wind);
        
        if nbr_regions_in_window >0         
            %fprintf("new window with %d regions in it \n", nbr_regions_in_window);

            % apply the filter 
            [ more_W_PL, indx_to_rm] = doublon_filter( MserRegions,local_indx_in_this_window, compress_percentage, same_size_percentage,min_dist,indx_to_rm,close_enough,bg, 'intra'); 

            %[W_PL , f_1_r_regions] = f1m6(regions,plot_f_1,bg,compress_percentage, same_size_percentage,min_dist,area_size_1,area_size_2,slide_by_raw, slide_by_col,min_range,max_range,close_enough );
            
            W_PL= [W_PL; more_W_PL];
            
            %W_PL = unique(W_PL, 'rows');
            [~, swpl] = size(W_PL);
            if swpl ~= 0
                [~,ia, ~] = unique(cell2mat(W_PL(:,2)));
                W_PL = W_PL(ia,:);
            end 
            [size_wpl, ~] = size(W_PL);
            for a=size_wpl : -1 : 1
                if ismember(cell2mat(W_PL(a,2)),indx_to_rm )
                    W_PL(a, :) = []; 
                end
            end
            
           [pl_size ,~]  = size(W_PL);
       %     fprintf('PL size: %d \n', pl_size);
            
            

%            [PL, more_index_to_rm] = filter1_in_window(PL,MserRegions,local_indx_in_this_window, compress_percentage, same_size_percentage,min_dist);
%            more_index_to_rm = filter1_in_window(PL,MserRegions,local_indx_in_this_window, compress_percentage, same_size_percentage,min_dist);
            
            %[first , ~] =size(more_index_to_rm);
                
           % indx_to_rm = unique([[indx_to_rm],[more_index_to_rm]]);

        end
        start_col= start_col + slide_by_col/2;
    end
    
    
    start_raw= start_raw + slide_by_raw/2;
end
        
        
%W_PL take only PL
PL =W_PL;
%[~,idu] = unique(PL(:,2));
%PL = PL(idu,:);
[~,ia, ~] = unique(cell2mat(PL(:,2)));
PL = PL(ia,1);

%PL = unique(PL,2);
%PL(:,2)= [];
%PL = unique(PL, 'rows');



[q,~] = size(PL);
%fprintf(' originally %d regions \n',MserRegions.Count); 
%fprintf(' PL size of keeped %d \n', q);
[~, tt] = size(indx_to_rm);
%fprintf(' %d  removed with filter_1 \n', tt);
%fprintf('difference : %d \n' , q-tt);

f_1_r_regions = MSERRegions(PL);





%PLOTTING
if plot_f_1
    %all_3_regr
    figure();imshow(bg); hold on;
    hp = impixelinfo;
    %figure; imshow(rgb); hold on;
    %plot(f_1_r_regions);
    %    plot(f_1_r_regions,'showPixelList',true);
    plot(f_1_r_regions,'showPixelList',false,'showEllipses',true);

end