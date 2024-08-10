function arr_sep = arrsepoptim(arr, gap, cirbnd)
%ARRSEPOPTIM Seperate the array for every "circular" distance between
% any two points are greater than gap. Make sure that the input range is
% [-cirbnd, cirbnd]. If not, the distance between two points may be
% negative. (cirbnd = +Inf by default). For cirbnd < 0, then assume that
% there is no circular bound, i.e. cirbnd = +Inf.
%
% The function simply makes use of the following MATLAB's built-in
% optimation toolbox:
% Solve the following optimization problem:
%   minimize    norm(v)^2
%   s.t.        A * (arr + v) >= gap * 1
%

if nargin < 3
	cirbnd = -1;
end

Q = length(arr);

if Q == 1
	arr_sep = arr;
elseif (cirbnd > 0) && (gap * Q > 2 * cirbnd)
	arr_sep = arr;
	warning("No feasible solution can be found.");
else
    %% check input
    if ~isvector(arr)
        error('Input must be a 1D array');
    end
    
    if ~iscolumn(arr)
        arr = arr.';
        flag_row = 1;
    else
        flag_row = 0;
    end

    %% check whether or not to compute optimization problem since it is time-consuming
    [arr, Ind] = sort(arr);

    if cirbnd > 0
        flag_cir = (2 * cirbnd - (arr(end) - arr(1))) < gap;
    else
        flag_cir = 0;
    end

    if any(diff(arr) < gap) || flag_cir
        prob = optimproblem("ObjectiveSense", "minimize");
	    v = optimvar('v', Q, 1);
    
	    % objective function
	    prob.Objective = norm(v)^2;
    
	    % constraints
	    A = full(spdiags([-ones(Q, 1), ones(Q, 1)], [0, 1], Q-1, Q));
    
        % 1. constraints between adjacent point
	    cons1 = A * (arr + v) >= gap;
	    prob.Constraints.cons1 = cons1;

        if cirbnd > 0
            % 2. constraints of the two end points
	        cons2 = arr(1) - arr(Q) + v(1) - v(Q) + 2 * cirbnd >= gap;
	        prob.Constraints.cons2 = cons2;

            % 3. boundary constraints
% 	        prob.Constraints.cons3 = arr(1) + v(1) >= -cirbnd;
% 	        prob.Constraints.cons4 = arr(Q) + v(Q) <= cirbnd;
        end
    
	    % solve the problem and unsort
	    options = optimoptions('lsqlin', 'Display', 'off');
	    [sol, ~, exitflag] = solve(prob, 'Solver', 'lsqlin', 'Options', options);
	    arr = arr + sol.v;
    
        if exitflag == -2
		    warning("No feasible solution found.");
        end

    end

    % if cirbnd > 0, deal with the circular boundary
    if cirbnd > 0
        for q = 1 : Q
            if arr(q) < -cirbnd
                arr(q) = arr(q) + 2 * cirbnd;
            end
            
            if arr(q) > cirbnd
                arr(q) = arr(q) - 2 * cirbnd;
            end
        end
    end

    % unsort
    arr_sep(Ind) = arr;

    if flag_row == 0
        arr_sep = arr_sep.';
    end

    % unsort
    arr_sep(Ind) = arr;

end

end

