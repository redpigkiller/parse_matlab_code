function arr = nan2zero(arr, val)
%NAN2ZERO Replace NaN in an array with val (0 by default)
%

if nargin < 2
    val = 0;
end

arr(isnan(arr)) = val;

end

