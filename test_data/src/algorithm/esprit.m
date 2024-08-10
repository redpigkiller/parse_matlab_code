function [f, val] = esprit(X, Q, args)
%ESPRIT This is the standard ESPRIT algorithm.
%
% X: received matrix
% Q: the number of estimated phase
% args:
%   sig_subspace:
%       SVD: square root approach
%       EVD: covariance approach
%   method:
%	    LS: (Least-Squares Solution)
%	    TLS: (Total-Least-Squares Optimization)
%

r = size(X, 1);
if Q > r;	error("can NOT solve! since Q > P");	end

%% performing ESPRIT

% 1. calculate U_hat
if strcmpi(args.sig_subspace, 'SVD')
    % i. compute SVD of T (square root approach):
    [U_hat, val] = svd(X);

elseif strcmpi(args.sig_subspace, 'EVD')
    % ii. compute EVD of T * T' (covariance approach):
    [U_hat, val] = eig(1/M * X * (X'));
	val = diag(val);
    [~, ind] = sort(val, 'descend');
	U_hat = U_hat(:, ind);

else
	error("wrong sig_subspace");
end

% 2. get U_hat_s
U_hat = U_hat(:, 1 : Q);

% 3. compute U_s1, U_s2
U_s1 = U_hat(1 : r - 1, :);
U_s2 = U_hat(2 : r, :);

% 4. choose method
if strcmpi(args.method, 'LS')
	% LS solution
	psi_LS = (U_s1'*U_s1)\U_s1'*U_s2;

	f = angle(eig(psi_LS)) / (2*pi);

elseif strcmpi(args.method, 'TLS')
	% TLS solution
	[~, ~, U] = svd([U_s1 U_s2]);
	U12 = U(1 : Q, Q+1 : end);
	U22 = U(Q+1 : end, Q+1 :end);

	psi_TLS = -U12/U22;
	
	f = angle(eig(psi_TLS)) / (2*pi);
else
	error("wrong method");
end

end