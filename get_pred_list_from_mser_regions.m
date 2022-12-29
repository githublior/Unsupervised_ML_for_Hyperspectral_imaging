function pred_list = get_pred_list_from_mser_regions(regions);
pred_list = [];
for region_index=1:regions.Count
    pred_list = [pred_list; {regions.PixelList(region_index)}];
end
