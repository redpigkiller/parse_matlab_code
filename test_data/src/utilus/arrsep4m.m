function new_arr = arrsep4m(arr, gap, cirbnd, tol)
%ARRSEP4M Seperate the array for every "circular" distance between
% any two points are greater than gap. Make sure that the input range is
% [-cirbnd, cirbnd]. If not, the distance between two points may be
% negative. (cirbnd = +Inf by default). For cirbnd < 0, then assume that
% there is no circular bound, i.e. cirbnd = +Inf.
%

%% check input arguments
if nargin < 3;  cirbnd  = -1;   end
if nargin < 4;  tol   = 1e-15;  end

assert(isvector(arr), 'Input must be a 1D array');
flag_row = ~iscolumn(arr);

%% main

% check for simple input
N = length(arr);
if N == 1
	new_arr = arr;
    return;
elseif (cirbnd > 0) && (gap * N > 2 * cirbnd)
	new_arr = arr;
	warning("No feasible solution can be found.");
    return;
end

% if not that simple
new_arr     = [arr(1)];
sort_arr    = [arr(1)];

for n = 2 : N
    new_arr     = insert_sort(new_arr, arr(n));
    sort_arr    = insert_sort(sort_arr, arr(n));

    % 1. linear part: find points which are too close
    d_arr   = diff(new_arr);
    Ind     = find(d_arr < gap);

    % 2. circular part: find if the two end points are too close
    if cirbnd > 0
        flag_cir = (2*cirbnd - (new_arr(n) - new_arr(1))) < gap;
    else
        flag_cir = 0;
    end

    while(~isempty(Ind) || flag_cir) 
        if ~isempty(Ind)
            % choose the first gap that is too close
            Ind = Ind(1);

            % the "distance" needed to be inserted
            d0 = gap - d_arr(Ind);

            % find contiguous left-hand side points
            count_l = 1;
            for idx = Ind - 1 : -1 : 1
                if (d_arr(idx) - gap) < tol
                    count_l = count_l + 1;
                else
                    break;
                end
            end
    
            % find contiguous right-hand side points
            count_r = 1;
            for idx = Ind + 1 : n-1
                if (d_arr(idx) - gap) < tol
                    count_r = count_r + 1;
                else
                    break;
                end
            end
    
            % calculate the "push distance"
            push_l = d0 * count_r / (count_r + count_l);
            push_r = d0 * count_l / (count_r + count_l);

            % deal with small value
            if push_l < eps(new_arr(Ind));  push_l = eps(new_arr(Ind)); end
            if push_r < eps(new_arr(Ind));  push_r = eps(new_arr(Ind)); end

            % push left
            for c = 1 : count_l
                new_arr(Ind+1-c) = new_arr(Ind+1-c) - push_l;
            end

            % push right
            for c = 1 : count_r
                new_arr(Ind+c) = new_arr(Ind+c) + push_r;
            end

        else
            
            % the "distance" needed to be inserted
            d1 = gap - (2*cirbnd - (new_arr(n) - new_arr(1)));
    
            % find contiguous left-hand side points
            count_l = 1;
            for idx = n-1 : -1 : 1
                if (d_arr(idx) - gap) < tol
                    count_l = count_l + 1;
                else
                    break;
                end
            end
    
            % find contiguous right-hand side points
            count_r = 1;
            for idx = 1 : n-1
                if (d_arr(idx) - gap) < tol
                    count_r = count_r + 1;
                else
                    break;
                end
            end
    
            % calculate the "push distance"
            push_l = d1 * count_r / (count_r + count_l);
            push_r = d1 * count_l / (count_r + count_l);

            % deal with small value
            if push_l < eps(new_arr(Ind));  push_l = eps(new_arr(Ind)); end
            if push_r < eps(new_arr(Ind));  push_r = eps(new_arr(Ind)); end
    
            % push left
            for c = 1 : count_l
                new_arr(n+1-c) = new_arr(n+1-c) - push_l;
            end
    
            % push right
            for c = 1 : count_r
                new_arr(c) = new_arr(c) + push_r;
            end
            
        end

        %% update
        % 1. linear part: find points which are too close
        d_arr = diff(new_arr);
        Ind = find(d_arr < gap);
    
        % 2. circular part: find if the two end points are too close
        if cirbnd > 0
            flag_cir = (2*cirbnd - (new_arr(n) - new_arr(1))) < gap;
        else
            flag_cir = 0;
        end
    end

end

% if cirbnd is positive, deal with the circular boundary
if cirbnd > 0
    for n = 1 : N
        if new_arr(n) < -cirbnd
            new_arr(n) = new_arr(n) + 2*cirbnd;
        end
        
        if new_arr(n) > cirbnd
            new_arr(n) = new_arr(n) - 2*cirbnd;
        end
    end
end

% re-order the new array
Ind     = arrcmp(arr, sort_arr);
new_arr = new_arr(Ind);

if flag_row
    new_arr = new_arr.';
end

end

function y = insert_sort(x, val)
N = length(x);

if N == 0
    y = val;
    return;
end    

% find the index of insertion
idx = 1;
while (idx <= N)
    if val < x(idx)
        break;
    end
    idx = idx + 1;
end

% insert the value to x
if idx == 1
    y = [val; x];
elseif idx == N + 1
    y = [x; val];
else
    y = [x(1:idx-1); val ; x(idx:end)];
end

end