function [f_2_regions] = filter_2_mser(wkn_img,MserRegions, void_in_ellipse_acceptance_rate,plot_f_2,bg,void_definition_rate,showPixelList,showEllipses,show_r_to_rm);

% percentage of miss to reject the element
indx_to_rm = [];
nbr_elt = MserRegions.Count;

for i=1:nbr_elt
    [pxl_grp,~] = size(MserRegions.PixelList(i));
    %fprintf('new group %d  \t there are %d pixels in this group \n', i,pxl_grp);

    miss_pixel=0;
    for j=1:pxl_grp
        %fprintf('new pixel %d \n',j);
        pxl_of_elt = MserRegions.PixelList(i);
        y=(pxl_of_elt(j,1));
        x=(pxl_of_elt(j,2));
        %wkn_img(x,y)
        if wkn_img(x,y)  < void_definition_rate
            miss_pixel = miss_pixel + 1 ;
        %else fprintf('pass');
        end
    end
    void_rate = double(miss_pixel)/pxl_grp;

    if void_rate > void_in_ellipse_acceptance_rate
        indx_to_rm =[indx_to_rm, i];
        
    end
    %size(indx_to_rm)
end

%removing them
PixelList = [];
pxl_to_rm = [];
for i=1:nbr_elt
    if  ~ismember(i, indx_to_rm)
        PixelList = [PixelList ;{MserRegions(i).PixelList}];
    else 
        pxl_to_rm = [pxl_to_rm ;{MserRegions(i).PixelList}];
    end
end
[~,nn]=size(indx_to_rm);
fprintf('%d objects removed with filter_2 \n', nn);
f_2_regions = MSERRegions(PixelList);


%PLOTTING

if show_r_to_rm
    region_to_rm = MSERRegions(pxl_to_rm);

    %all_3_regr
    figure('Name','removed area ');imshow(bg); hold on;
    hp = impixelinfo;
    %figure; imshow(all_3_regr); hold on;
    %plot(f_2_regions);
    plot(region_to_rm,'showPixelList', showPixelList, 'showEllipses',showEllipses);
end
if plot_f_2
    %all_3_regr
    figure('Name','filter_2 ');imshow(bg); hold on;
    hp = impixelinfo;
    %figure; imshow(all_3_regr); hold on;
    %plot(f_2_regions);
    plot(f_2_regions,'showPixelList', showPixelList, 'showEllipses',showEllipses);
end
