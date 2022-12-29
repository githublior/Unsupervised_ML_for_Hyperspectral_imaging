function [mean_h_for_this_region, mean_e_for_this_region] =get_hne_mean_of_mser_regions_pixellist(regions_pixels, H, E);
    [s,~] = size(regions_pixels);
    indv_h_val =0;
    indv_e_val =0;
    for coord=1:s
        pxl_val_in_h =uint32(H(regions_pixels(coord,2), regions_pixels(coord,1)));
        indv_h_val = indv_h_val +pxl_val_in_h;
        pxl_val_in_e = uint32(E(regions_pixels(coord,2), regions_pixels(coord,1)));
        indv_e_val = indv_e_val + pxl_val_in_e;
    end
    mean_h_for_this_region = indv_h_val/s;
    mean_e_for_this_region = indv_e_val/s;
end