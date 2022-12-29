function pixel_nbr_in_cell = get_cell_size_from_user(rgb, cell_type);

    if cell_type == 'large'
        fprintf('label a large cell in the tissue.\n');
    else
        fprintf('label a small cell in the tissue.\n');
    end
    imshow(rgb);
    large_cell = drawellipse();
    %{
y= 'n';
    while y ~= 'y'
        imshow(rgb);
        large_cell = drawellipse();
        y = input('Press y when done.. \n','s');
    end
    %}
    large_cell_mask = createMask(large_cell);
    pixel_nbr_in_cell = sum(sum(large_cell_mask));


end