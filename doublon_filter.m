%more_index_to_rm
%close_enough = 20;

function  [ W_PL,more_index_to_rm] = doublon_filter( regions, local_indx_in_this_window, compress_percentage, same_size_percentage,min_dist, indx_to_remove, close_enough, bg, filter_type);
% define min_dist according to hafifot.
dynamic_local_index =local_indx_in_this_window;
%apply classic filter_1 to the local mserRegions.
W_PL = [];
more_index_to_rm = indx_to_remove;

[s_loc, ~] = size(dynamic_local_index);
i =1;
j=i+1;

while i<=s_loc
    if i==111
        ke=1;
    end
    if j > s_loc
        [s_d_l_i, ~] = size(dynamic_local_index);
        if s_d_l_i == 0 || ~ismember(dynamic_local_index(i) , more_index_to_rm)
            W_PL = [W_PL; {regions(dynamic_local_index(i)).PixelList} , [dynamic_local_index(i)] ];
        end
        i = i+1;
        j = i+1;
        continue
    end
    
    
    %get loc x y of mser i         
    loc_xandy_i = regions(dynamic_local_index(i)).Location;
    x_l_i =double(loc_xandy_i(1,1));
    y_l_i =double(loc_xandy_i(1,2));

    %get loc x y of mser j

    loc_xandy_j = regions(dynamic_local_index(j)).Location;
    x_l_j =double(loc_xandy_j(1,1));
    y_l_j =double(loc_xandy_j(1,2));
    

    if  ismember(dynamic_local_index(i),more_index_to_rm)== 1 
        i = i+1;
        j = i+1;
        continue
    end    
    if ismember(dynamic_local_index(j),more_index_to_rm)== 1
        j = j+1;
        continue
    end 
    
    
    
      
    x_low_bound = 150;
    x_high_bound = 160;
    y_low_bound = 165;
    y_high_bound = 175;
    if x_l_i > x_low_bound && x_l_j > x_low_bound && y_l_i> y_low_bound &&  y_l_j> y_low_bound
        if x_l_i < x_high_bound && x_l_j  < x_high_bound && y_l_i < y_high_bound &&  y_l_j < y_high_bound
            tetu=1;
        end
    end
    
    
    
    %case the center of both regions are close enough to be the same area   
    if abs(x_l_i - x_l_j) < close_enough && abs(y_l_i - y_l_j) < close_enough 
        [size_i,~]  = size(regions(dynamic_local_index(i)).PixelList);
        [size_j,~]  = size(regions(dynamic_local_index(j)).PixelList);
        
        if filter_type == 'intra'  %then select the largest area (nv)
            if max([size_i , size_j]) == size_i
                more_index_to_rm = [more_index_to_rm, dynamic_local_index(j)];
                j  = j+1;
            else
                more_index_to_rm = [more_index_to_rm, dynamic_local_index(i)];
                i = i+1;
                j = i+1;
            end
        else   %then select the smallest area (nv)
            if min([size_i , size_j]) == size_i
                more_index_to_rm = [more_index_to_rm, dynamic_local_index(j)];
                j  = j+1;
            else
                more_index_to_rm = [more_index_to_rm, dynamic_local_index(i)];
                i = i+1;
                j = i+1;
            end
            
        end  

        
    %case i and j not same area    
    else
        j = j+1;
             
    end
    
    [s_loc, ~] = size(dynamic_local_index);
    
end

%reduce W_PL to unique
[~, swpl] = size(W_PL);
if swpl ~= 0
    [~,ia, ~] = unique(cell2mat(W_PL(:,2)));
    W_PL = W_PL(ia,:);
end

end
%old filter 1 com
%{
for i=1:s_loc
    for j=i+1:s_loc
        if i==j 
            continue;

    full =  intersect(regions(local_indx_in_this_window(i)).PixelList,regions(local_indx_in_this_window(j)).PixelList,'rows');
    
        [f1,~] = size(full);
        


        [p,~] = size(W_PL);
        
        %case first line of W_PL
        if p==0
            % case ~i and ~j are different
            if f1 ==0
                W_PL = [W_PL; {regions(local_indx_in_this_window(i)).PixelList} , [local_indx_in_this_window(i)] ];
                W_PL = [W_PL; {regions(local_indx_in_this_window(j)).PixelList} , [local_indx_in_this_window(j)] ];
            
            % case ~i and ~j are same cell
            else
                merge = cat(1,regions(local_indx_in_this_window(i)).PixelList,regions(local_indx_in_this_window(j)).PixelList);
                W_PL = [W_PL; {unique(merge,'rows')} , [local_indx_in_this_window(i),local_indx_in_this_window(j)] ];
            end
        
        else
            
            

            i_indx =-1;
            j_indx =-1;
            i_touch =false;
            j_touch=false;

            for ch=1:p
                %check if  this row contain index ~i in the array (2nd col)
                idxi = find(W_PL{ch , 2}==local_indx_in_this_window(i));
                [~ ,i2] = size(idxi);

                if i2>0
                    i_indx = ch;
                    i_touch=true;
                end

                idxj=find(W_PL{ch , 2}==local_indx_in_this_window(j) );
                [~ ,j2]= size(idxj);

                if j2>0
                    j_indx = ch;
                    j_touch =true;
                end

                if i_touch && j_touch
                    break
                end

            end

            
            
            %if they have pixels in common
            if f1  >= 1 
                %list that indx PL with indx ~i , ~j such that if we have a new
                %merge including ~i or ~j , we merge with the previous instead
                %of a whole new line.





                %CASE none of them are in window_PL:create 1 new line for 2
                if ~i_touch && ~j_touch
                    %update more_index_to_rm
                    more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(i)];
                    more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(j)];
                    %create simple new line
                    merge=cat(1,regions(local_indx_in_this_window(i)).PixelList,regions(local_indx_in_this_window(j)).PixelList);
                    W_PL= [W_PL ; {unique(merge,'rows')} , [local_indx_in_this_window(i), local_indx_in_this_window(j)] ];
                end

                %CASE both ~i and ~j are found in window_PL . 2 cases: in same
                %or different lines.
                if i_touch && j_touch
                    %case same line in window_PL: dont add anything
                    if i_indx == j_indx
                        continue;

                    %case different lines     
                    else
                        %merge les 2 lignes  de window_PL ds la 1ere et rm la deuxieme

                        %merge the pxllist in ~i
                        merge_in_w_pl=cat(1,W_PL{i_indx,1},W_PL{j_indx,1});
                        W_PL(i_indx, 1)  = {unique(merge_in_w_pl,'rows')};


                        %merge indx array in ~i 
                        W_PL{i_indx,2} = unique([[W_PL{i_indx,2}] , [W_PL{j_indx,2}]]);


                        %rm j line since mergin with i 
                        W_PL(j_indx,:) = [];

                        %update more_index_to_rm
                        more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(i)];
                        more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(j)];
                        continue


                    end

                end




                %CASE  seen ~i but not ~j : add j to line i_indx
                if i_touch && ~j_touch         


                    %change the pxllist of i to add j pixllist
                    merge_in_w_pl=cat(1,W_PL{i_indx,1},regions(local_indx_in_this_window(j)).PixelList);
                    W_PL(i_indx, 1) = {unique(merge_in_w_pl,'rows')};

                    %update the indx array i with  j 
                    W_PL{i_indx,2} = unique([[W_PL{i_indx,2}], [local_indx_in_this_window(j)]]);


                    %update more_indx_to_rm
                    more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(j)];

                    continue
                end



                %CASE OPPOSITE : not~i but ~j :  add ~i to line j_idx
                if ~i_touch && j_touch         

                    %change the pxllist of j to add i pixllist
                    merge_in_w_pl=cat(1,W_PL{j_indx,1},regions(local_indx_in_this_window(i)).PixelList);
                    W_PL(j_indx, 1) = {unique(merge_in_w_pl,'rows')};

                    %update the indx array j with  i 
                    W_PL{j_indx,2} = unique([[W_PL{j_indx,2}], [local_indx_in_this_window(i)]]);


                    %update more_indx_to_rm
                    more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(j)];

                    continue
                end            



            else
            %CASE f1=0 ie the 2 regions doesnt have pxl in common.    
            %now check wether ~i and ~j are already in W_PL:
  

                %Case there are not the same area (no CROSSING)

                if ~i_touch && ~j_touch
                    





                end


                % CASE not ~i and not ~j add ~i and ~J
                if   ~i_touch && ~j_touch   
                    %update more_indx_to_rm
                    more_index_to_rm =[more_index_to_rm, loca
                        l_indx_in_this_window(i)];
                    more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(j)];

                    %add new line in w_pl
                    W_PL=[W_PL; {regions(local_indx_in_this_window(i)).PixelList}, [local_indx_in_this_window(i)]];
                    W_PL=[W_PL; {regions(local_indx_in_this_window(j)).PixelList}, [local_indx_in_this_window(j)]];

                end
                % CASE ~I AND NOT ~J : ADD ~J ONLY 
                if i_touch && ~j_touch
                    %update more_indx_to_rm
                    more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(j)];
                    %add new line in w_pl
                    W_PL=[W_PL; {regions(local_indx_in_this_window(j)).PixelList}, [local_indx_in_this_window(j)]];
                end
                % CASE NOT ~I AND  ~J : ADD ~I ONLY 
                if ~i_touch && j_touch
                    %update more_indx_to_rm
                    more_index_to_rm =[more_index_to_rm, local_indx_in_this_window(i)];
                    %add new line in w_pl
                    W_PL=[W_PL; {regions(local_indx_in_this_window(i)).PixelList}, [local_indx_in_this_window(i)]];
                end


            end
            more_index_to_rm =unique(more_index_to_rm);

        end

    end

end
%}

%printing helper
%{
cri =   MSERRegions({regions(dynamic_local_index(i)).PixelList});
crj =   MSERRegions({regions(dynamic_local_index(j)).PixelList});imshow(bg); hold on;
hp = impixelinfo;
plot(cri,'showPixelList',false,'showEllipses',true);
plot(crj,'showPixelList',false,'showEllipses',true);

%}