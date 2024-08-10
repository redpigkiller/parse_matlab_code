function Xcir = circulant(x, Lch)
%CIRCULANT Generate circulant matrix Xcir from vector x.
%
% x:	input vector
% Lch:	number of circulant shift
% generate N x Lch matrix

	N = size(x, 1);

	Xcir = zeros(N, Lch, 'like', x);
	for i_Lch = 1 : Lch
		Xcir(:, i_Lch) = circshift(x, i_Lch-1);
	end
end