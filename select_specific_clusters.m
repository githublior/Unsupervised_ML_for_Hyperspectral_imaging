function [bg] = select_specific_clusters( bg, working_img, cell_area);
    %{
    %for true cell analysis
    blank: .22, .56, .67 little .33
    cyto : .78 1 , .89
    kernel/membrane : .44
    between kernel-cyto :0 
    %}

    % miss between kernel cytoplasm : 57 - 85
    [a,b] = size(working_img);
    [~,c] = size(cell_area);

    for i = 1:a
        for j= 1:b
            bol =0;
            for k= 1:c
                if uint8(working_img(i,j)) == cell_area(k)
                    bol = 1;
                end
            end
            if bol ==1
                bg(i,j)=1;
            end
            %else
             %   binay(i,j) = 0; 
                %binay(i,j) = mat2gray(kernel_out(i,j));
                %binay(i,j) = max(mat2gray(kernel_out(i,j))-.5,0);
        end
    end
