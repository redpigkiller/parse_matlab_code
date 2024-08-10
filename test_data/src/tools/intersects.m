function [tgt, unmatched] = intersects(src, varargin)
%INTERSECTS Compares multiple input arrays, return their intersect
%
%   See also INTERSECT.

unmatched = 0;

if nargin <= 1
	tgt = src;

else
    if length(src) ~= length(varargin{1});  unmatched = 1;  end
	tgt = intersect(src, varargin{1});

	for ii = 2 : nargin - 1
        if length(tgt) ~= length(varargin{ii});  unmatched = 1;  end
		tgt = intersect(tgt, varargin{ii});
	end
end

end

