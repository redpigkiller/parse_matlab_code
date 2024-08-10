function [arr_hat, indices, cost] = arrpairbruteforce(arr_ref, arr_hat, cirbnd)
%ARRPAIRBRUTEFORCE Compute optimal assignment by iterating through all
% possible assignment.
% The first input arr_ref is the reference array that arr_hat to pair to.
% The value of the ouput is the re-order of the second input arr_hat.
% Make sure the input range is [-cirbnd, cirbnd], otherwise, the distance
% between some two points may be negative. (cirbnd = +Inf by default)
%

%% check
if ~isvector(arr_ref) || ~isvector(arr_hat)
    error('Inputs must be a 1D array');
end

if ~iscolumn(arr_hat)
    arr_hat = arr_hat.';
    flag_row = 1;
else
    flag_row = 0;
end

if nargin < 3
    cirbnd = -1;
end

%% main
Q = length(arr_ref);

% distance matrix
dist = zeros(Q, Q, 'like', arr_ref);

if cirbnd > 0
    for i_Q = 1 : Q
	    for ii_Q = 1 : Q
            d = abs(arr_ref(i_Q) - arr_hat(ii_Q));
		    dist(i_Q, ii_Q) = min(d, 2 * cirbnd - d);
	    end
    end
else
    for i_Q = 1 : Q
	    for ii_Q = 1 : Q
		    dist(i_Q, ii_Q) = abs(arr_ref(i_Q) - arr_hat(ii_Q));
	    end
    end
end

% brute-force
Ind = perms(1:Q);
indices = size(Ind, 1);

cost = inf;
min_idx = 0;
for p = 1 : indices
	tot_cost = 0;
	for q = 1 : Q
		tot_cost = tot_cost + dist(q, Ind(p, q));
	end
	if tot_cost < cost
		cost = tot_cost;
		min_idx = p;
	end
end

% re-order
indices = Ind(min_idx, :).';
arr_hat = arr_hat(indices);

if flag_row == 1
    arr_hat = arr_hat.';
    indices = indices.';
end

end