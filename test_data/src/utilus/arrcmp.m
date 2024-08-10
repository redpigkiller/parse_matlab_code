function indices = arrcmp(sim_arr, exist_arr, del_zeros)
%ARRCMP Compares two input arrays, return the indices corresponding to the
%exist_arr
%
%   Example1: sim_arr = [1, 3, 5] and exist_arr = [1, 2, 3, 4, 5]
%       indices = [1, 3, 5]
%   Example2: sim_arr = [1, 2, 3, 4, 5] and exist_arr = [1, 3, 5]
%       indices = [1, 0, 2, 0, 3]
%

if nargin < 3
    del_zeros = 1;
end

n1 = length(sim_arr);
n2 = length(exist_arr);

indices = zeros(n1, 1);
for ii = 1 : n1
	for jj = 1 : n2
		if sim_arr(ii) == exist_arr(jj)
			indices(ii) = jj;
			break;
		end
	end
end

if del_zeros ~= 0
	indices(indices == 0) = [];
end

end
