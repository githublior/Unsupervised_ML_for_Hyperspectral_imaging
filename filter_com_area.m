%more_index_to_rm
function  [ W_PL,more_index_to_rm] = filter_com_area( regions, local_indx_in_this_window, compress_percentage, same_size_percentage,min_dist, indx_to_remove);
% define min_dist according to hafifot.


%apply classic filter_1 to the local mserRegions.
W_PL = [];
more_index_to_rm = indx_to_remove;

[s_loc, ~] = size(local_indx_in_this_window);

for i=1:s_loc
    %print run time
    %fprintf('%d / %d \n',i, s_loc);
    
    %if ismember(local_indx_in_this_window(i),[436,437]) 
    %   nothing=0;
    %   fprintf('%d %d \n', i ,local_indx_in_this_window(i) );
   % end
   
   % if ismember(local_indx_in_this_window(i), more_index_to_rm) 
   %     %fprintf('for i = %d , loc-ondx_to_rm(i)= %d ,and alreday in more_indx_to_rm  \n' , i ,local_indx_in_this_window(i) );
   %     continue
  %  end


    for j=i+1:s_loc
        %if ismember(local_indx_in_this_window(j),[436,437]) 
          %  nothing = 0 ;
         %   fprintf('%d %d \n', j ,local_indx_in_this_window(j) );
        %end
        
        %if ismember(local_indx_in_this_window(j), more_index_to_rm)
        %    %fprintf('for j = %d , loc-ondx_to_rm(j)= %d ,and alreday in more_indx_to_rm  \n' , j ,local_indx_in_this_window(j) );
        %    continue
        %end
        
        
        %fprintf('i= %d  ,local_indx_in_this_window(i)= %d , j=%d ,local_indx_in_this_window(j)=%d both are not in more_indx_to_rm \n  ', i,local_indx_in_this_window(i) , j , local_indx_in_this_window(j) );
        
        
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



