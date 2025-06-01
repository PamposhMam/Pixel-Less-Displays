function [ordered,leftover] = ordered(A)      
    points = length(A(:,1));   
    
    k = 2;
    distances = [];
    ordered = [A(1,:)];
    
    x_orig = A(1,1); y_orig = A(1,2);   %store initial point
    x = x_orig; y = y_orig;             %initialize for algorithm
    num_search = size(A(:,1));
    while length(A(:,1)) > 0
        k1 = k;                         %store index of initial point
    
        while abs(x - A(k,1)) <= 2e5      %search until points are far from initial point in x
            dist = sqrt( (x - A(k,1))^2 + (y - A(k,2))^2 );     %compute distance from initial point
            distances = [distances; k, dist];                   %array of distances and point index
            k = k+1;                                            %search the next point
            if k > num_search                                   %break if at the end of the array
                break;
            end
        end
    
        k = k1-1;                       %search from initial point in other direction
        if k > 0                        %if not at the start of the array
            while abs(x - A(k,1)) <= 2e5 && k > 0                 %search until points are far from initial point in x && not at the start of array
                dist = sqrt( (x - A(k,1))^2 + (y - A(k,2))^2 ); %compute distance from initial point
                distances = [distances; k, dist];               %array of distances and point index
                k = k-1;                                        %search the next point
                if k <= 0                                       %break if at the end of the array
                    break;
                end
            end
        end
    
        [min_dist,index] = min(distances(:,2));                 %out of all the nearby points, which is closest?
        if min_dist < 10                                        %if not at the end of the curve
            index = distances(index,1);                         
            ordered = [ordered; A(index,:)];                    %add the nearest point to the end of the ordered array
            x = A(index,1);                                     %save this point as the new initial point for next search
            y = A(index,2);
            A = [A(1:index-1,:);A(index+1:num_search,:)];       %remove the point from the array of to-be-sorted points
            k = index;                                          %
            distances = [];                                     %clear the array of distances computed
            num_search = length(A(:,1));    %find max number of points to search
        else
            break;
        end
    
        if k > num_search           %if the end of the array is reached
            break;                  %end the loop
        end
    end 
    
    % check ordered(1) vs A(1) and A(n) and ordered(n) vs A(1) and A(n)
    if length(A) > 0    % if there are still more points to order
    
        endpoints = [ordered(1,:); ...
            ordered(length(ordered),:); ...
            A(1,:);...
            A(length(A(:,1)),:)];
    
        check_endpoints = ...       %[1,3;1,4;2,3;2,4]
        [sqrt( (endpoints(1,1) - endpoints(3,1))^2 + (endpoints(1,2) - endpoints(3,2))^2 );...  find distance between the first element of the ordered vector and first element of remaining vector
        sqrt( (endpoints(1,1) - endpoints(4,1))^2 + (endpoints(1,2) - endpoints(4,2))^2 );...   
        sqrt( (endpoints(2,1) - endpoints(3,1))^2 + (endpoints(2,2) - endpoints(3,2))^2 );...
        sqrt( (endpoints(2,1) - endpoints(4,1))^2 + (endpoints(2,2) - endpoints(4,2))^2 )];
    
        [min_check_endpoints,idx] = min(check_endpoints);
        switch idx
            case 1                                          %match first element of sorted vector and first element of remaining vector
                k = 1;                                      %start searching first element of remaining points
                x = endpoints(1,1); y = endpoints(1,2);     %compare distance to first element of sorted vector
            case 2                                          %match first element of sorted vector and last element of remaining vector
                k = length(A(:,1));                         %start searching last element of remaining points
                x = endpoints(1,1); y = endpoints(1,2);     %compare distance to first element of sorted vector
            case 3                                          %match last element of sorted vector and first element of remaining vector
                k = 1;                                      %start searching first element of remaining points
                x = endpoints(2,1); y = endpoints(2,2);     %compare distance to last element of sorted vector
            case 4                                          %match last element of sorted vector and last element of remaining vector
                k = length(A(:,1));                         %start searching last element of remaining points
                x = endpoints(2,1); y = endpoints(2,2);     %compare distance to last element of sorted vector
        end
        num_search = length(A(:,1));                        %find max number of points to search
    
        while length(A) > 0                                 %until all points have been sorted
            k1 = k;                                         %store index of initial point
    
            while abs(x - A(k,1)) <= 2e5                      %search until points are far from initial point in x
                dist = sqrt( (x - A(k,1))^2 + (y - A(k,2))^2 );     %compute distance from initial point
                distances = [distances; k, dist];                  %array of distances and point index
                k = k+1;                                            %search the next point
                if k > num_search                                   %break if at the end of the array
                    break;
                end
            end
    
            k = k1-1;                                               %search from initial point in other direction
            if k > 0                                                %if not at the start of the array
                while abs(x - A(k,1)) <= 2e5 && k > 0                 %search until points are far from initial point in x && not at the start of array
                    dist = sqrt( (x - A(k,1))^2 + (y - A(k,2))^2 ); %compute distance from initial point
                    distances = [distances; k, dist];               %array of distances and point index
                    k = k-1;                                        %search the next point
                    if k <= 0                                       %break if at the end of the array
                        break;  
                    end
                end
            end
    
            [min_dist,index] = min(distances(:,2));                 %out of all the nearby points, which is closest?
            if min_dist < 100                                        %if not at the end of the curve
                index = distances(index,1);
                if idx <= 2
                    ordered = [A(index,:);ordered];                 %add the nearest point to the start of the ordered array
                elseif idx > 2
                    ordered = [ordered; A(index,:)];                %add the nearest point to the end of the ordered array
                end
                k = index;
                x = A(index,1);                                     %save this point as the new initial point for next search
                y = A(index,2);
                A = [A(1:index-1,:);A(index+1:num_search,:)];       %remove the point from the array of to-be-sorted points
    
                distances = [];                                     %clear the array of distances computed
                num_search = length(A(:,1));    %find max number of points to search
                if k > num_search               %make sure not to go out of array bounds
                    k = k-1;
                end
            else
                break;
            end
        end 
    end
    
    if ordered(1,1) > ordered(length(ordered(:,1)),1)       %first element of ordered vector
        ordered = flip(ordered);                            %is forced to be the most negative x-coord
    
    elseif ordered(1,1) == ordered(length(ordered(:,1)),1)  %if endpoints are the same x-coord
        if ordered(1,2) > ordered(length(ordered(:,2)),2)   %first element of ordered vector
            ordered = flip(ordered);                        %is forced to be the most negative y-coord
        end
    end

    leftover = A;
end