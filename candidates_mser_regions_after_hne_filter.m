function [candidates_regions_pxlist] = candidates_mser_regions_after_hne_filter(full_H, full_E,msers_candidates_regions, H_E_quantile, selected_area);


    H_acceptance_treshold_value = get_hne_acceptance_treshold( full_H, H_E_quantile);
    E_acceptance_treshold_value = get_hne_acceptance_treshold( full_E, H_E_quantile);

    [~ , s2_sa] = size(selected_area);
    if s2_sa == 0
        H = full_H;
        E = full_E;
    else
        H = imcrop(full_H, selected_area);
        E = imcrop(full_E, selected_area);
    end

    candidates_regions_pxlist = [];

    for region_index=1:msers_candidates_regions.Count
        regions_pixels = msers_candidates_regions.PixelList(region_index);
        [mean_h_for_this_region, mean_e_for_this_region] =get_hne_mean_of_mser_regions_pixellist(regions_pixels, H, E);
        if mean_h_for_this_region>= H_acceptance_treshold_value || mean_e_for_this_region>= E_acceptance_treshold_value
            candidates_regions_pxlist = [candidates_regions_pxlist;{msers_candidates_regions(region_index).PixelList}];
        end
    end
end