function []= manual_cell_labelization(tissue_name, img, crop_area, crop_area_name, H, E);

    x=0;
    b=false;
    

%qwert = insertShape(background_full,'filledrectangle', crop_area_1, 'Color', 'blue','Opacity',.3);
%qwert = insertShape(qwert,'filledrectangle', crop_area_2, 'Color', 'blue','Opacity',.3);
%qwert = insertShape(qwert,'filledrectangle', crop_area_3, 'Color', 'blue','Opacity',.3);
%qwert = insertShape(qwert,'filledrectangle', crop_area_4, 'Color', 'blue','Opacity',.3);

    while x<100 && b==false
        imtool(H);
        imtool(E);
        imshow(img);
        roi = drawellipse(); % select manually gt

        
        
        x = x+1; % put DEBUG here: if labelization over: write b=true in Matlab Terminal
        bw = createMask(roi);
        rd = get_radius_of_ellipse(roi);
        img = insertShape(img,'FilledCircle', [roi.Center(1) roi.Center(2) rd], 'Color', 'white','Opacity',.6);
        H = insertShape(H,'FilledCircle', [roi.Center(1) roi.Center(2) rd], 'Color', 'blue','Opacity',.3);
        E = insertShape(E,'FilledCircle', [roi.Center(1) roi.Center(2) rd], 'Color', 'blue','Opacity',.3);
        save('/Users/lior/Desktop/Image & Analysis /cells_label/'+ string(tissue_name) +'/'+ string(crop_area_name) + '/cell_' + string(num2str(x)) + '.mat', 'bw');
    end
end