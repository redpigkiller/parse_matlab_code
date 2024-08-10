function digits = baseN(num, base, num_digits)
% BASEN Convert a decimal number to a base-N number with a specified number of digits
% num:  the decimal number to convert
% base: the desired base
% num_digits (optional): the number of digits in the resulting base-N number
%   if not provided, calculates it automatically

if nargin < 3
    % if num_digits not provided, calculate it automatically
    num_digits = ceil(log(num+1) / log(base));
end

digits = zeros(1, num_digits);
for idx = num_digits : -1 : 1
    digits(idx) = mod(num, base);
    num = floor(num / base);
end

end

