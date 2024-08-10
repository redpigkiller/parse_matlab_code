function dist = circdist(arr1, arr2, circ_bnd, ignore_check)
%CIRCDIST Return the circular distance of arr1 and arr2.
%

if nargin < 4
    ignore_check = 0;
end

if ignore_check == 0
    assert(all(abs(arr1) <= circ_bnd) && all(abs(arr2) <= circ_bnd), ...
        "Array elements not in the range of [-cir_bnd, cir_bnd]");
end

abs_dist = abs(arr1 - arr2);
dist = min(abs_dist, 2 * circ_bnd - abs_dist);

end

