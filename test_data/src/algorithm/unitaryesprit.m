function [f, r_test] = unitaryesprit(X, Q, args)
%UNITARYESPRIT This is the unitary ESPRIT algorithm.
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

P = size(X, 1);
if Q > P;	error("can NOT solve! since Q > P");	end

%% 1. Real-Valued Subspace Estimation
% a. generate unitary left II-real matrix
if mod(P, 2) == 0
    QM = 1/sqrt(2) * [
        eye(P/2, 'like', X),         1j * eye(P/2, 'like', X)
        flip(eye(P/2, 'like', X)),   -1j * flip(eye(P/2, 'like', X))
        ];
    Qm = 1/sqrt(2) * [
        eye(P/2-1, 'like', X),       zeros(P/2-1, 1, 'like', X), 1j * eye(P/2-1, 'like', X)
        zeros(1, P/2-1, 'like', X),  sqrt(2),                    zeros(1, P/2-1, 'like', X)
        flip(eye(P/2-1, 'like', X)), zeros(P/2-1, 1, 'like', X), -1j * flip(eye(P/2-1, 'like', X))
        ];
else
    QM = 1/sqrt(2) * [
        eye((P-1)/2, 'like', X),       zeros((P-1)/2, 1, 'like', X), 1j * eye((P-1)/2, 'like', X)
        zeros(1, (P-1)/2, 'like', X),  sqrt(2),                      zeros(1, (P-1)/2, 'like', X)
        flip(eye((P-1)/2, 'like', X)), zeros((P-1)/2, 1, 'like', X), -1j * flip(eye((P-1)/2, 'like', X))
        ];
    Qm = 1/sqrt(2) * [
        eye((P-1)/2, 'like', X),        1j * eye((P-1)/2, 'like', X)
        flip(eye((P-1)/2, 'like', X)),  -1j * flip(eye((P-1)/2, 'like', X))
        ];
end

% b. the unitary transformation
%   (i)     compute T = Q_M' * [ X, II_M * conj(X) * II_N ] * Q_2N
T = 1/sqrt(2) * QM' * [X + flip(conj(X), 1), 1j * (X - flip(conj(X), 1))];

%   (ii)    T should be real
T = real(T);

% c. compute the signal subspace
%   (i)     choose method
if strcmpi(args.sig_subspace, 'SVD')
    % compute SVD of T (square-root approach):
    [Es, ~] = svd(T);

elseif strcmpi(args.sig_subspace, 'EVD')
    % compute EVD of T * T' (covariance approach):
    [Es, eigvar] = eig(T * T');
	eigvar = diag(eigvar);
    [~, ind] = sort(eigvar, 'descend');
	Es = Es(:, ind);

else
	error("wrong sig_subspace");
end

%   (ii)    select the Q dominated left singular/eigen vectors
Es = Es(:, 1:Q);

%% 2. Real-Valued Invariance Equation
% solve K1 * Es * R ~= K2 * Es
% a. preparation
J2 = [zeros(P-1, 1, 'like', X), eye(P-1, 'like', X)];
Qm = Qm' * J2 * QM;
K1 = 2 * real(Qm) * Es;
K2 = 2 * imag(Qm) * Es;

% choose method:
if strcmpi(args.method, 'LS')
    % a. least square (LS) solution:
	R = (K1'*K1)\K1'*K2;

elseif strcmpi(args.method, 'TLS')
    % b. total least square (TLS) solution:
	[~, ~, U] = svd([K1 K2]);
	U12 = U(1 : Q, Q+1 : end);
	U22 = U(Q+1 : end, Q+1 :end);

	R = -U12/U22;
	
else
	error("wrong method");
end

%% DOA Estimation
w = eig(R);

%% Reliability Test
if any(imag(w) ~= 0)
    % failed
    r_test = 0;
else
    % passed
    r_test = 1;
end
w = real(w);

%% Angle Estimation
f = atan(w) / pi;

end