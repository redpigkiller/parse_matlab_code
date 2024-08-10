function [fu, fv] = unitaryesprit_2d(X, Mx, Q, args)
%UNITARYESPRIT_2D This is the 2D unitary ESPRIT algorithm.
%
% The function assumes that the input is from an UPA (M = Mx * My).
% The input should be arranged in the following form:
% Denote the output from the (i, j) antenna y_{i,j} (Nx1) with N snapshots.
% Then y_{i,j} at time t are stacked in a column vector x(t) (M x 1).
%
% X:    [x(0) x(1) ... x(N-1)] is a M x N matrix
% Mx:   the number of antennas along the first dimension
% Q:    the number of estimated phase
% args:
%   sig_subspace:
%       SVD: square root approach
%       EVD: covariance approach
%   method:
%	    LS: (Least-Squares Solution)
%	    TLS: (Total-Least-Squares Optimization)
%

M = size(X, 1);
assert(mod(M, Mx) == 0, "invalid array manifold!");
My = M / Mx;

%% 1. Real-Valued Subspace Estimation
% a. generate unitary left II-real matrix
QM = pi_real(M);

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
    
%   (ii)    select the Q dominat left singular/eigen vectors
Es = Es(:, 1:Q);

%% 2. Real-Valued Invariance Equation
% solve Ku1 * Es * Ru ~= Ku2 * Es
% solve Kv1 * Es * Rv ~= Kv2 * Es
% a. preparation
Qmx = pi_real((Mx-1) * My);
Qmy = pi_real(Mx * (My-1));

Ju2 = [zeros(Mx-1, 1, 'like', X), eye(Mx-1, 'like', X)];
Jv2 = [zeros(My-1, 1, 'like', X), eye(My-1, 'like', X)];
Ju2 = kron(eye(My, 'like', X), Ju2);
Jv2 = kron(Jv2, eye(Mx, 'like', X));

Qmx = Qmx' * Ju2 * QM;
Qmy = Qmy' * Jv2 * QM;

Ku1 = 2 * real(Qmx) * Es;
Ku2 = 2 * imag(Qmx) * Es;
Kv1 = 2 * real(Qmy) * Es;
Kv2 = 2 * imag(Qmy) * Es;

% choose method:
if strcmpi(args.method, 'LS')
    % a. least square (LS) solution:
	Ru = (Ku1'*Ku1)\Ku1'*Ku2;
	Rv = (Kv1'*Kv1)\Kv1'*Kv2;

elseif strcmpi(args.method, 'TLS')
    % b. total least square (TLS) solution:
	[~, ~, U] = svd([Ku1 Ku2]);
	U12 = U(1 : Q, Q+1 : end);
	U22 = U(Q+1 : end, Q+1 :end);

	Ru = -U12/U22;

    [~, ~, U] = svd([Kv1 Kv2]);
	U12 = U(1 : Q, Q+1 : end);
	U22 = U(Q+1 : end, Q+1 :end);

	Rv = -U12/U22;
	
else
	error("wrong method");
end

%% DOA Estimation
% Automatic Pairing of Spatial Frequency Estimates
w = eig(Ru + 1j * Rv);

%% Angle Estimation
fu = atan(real(w)) / pi;
fv = atan(imag(w)) / pi;

end

function Qm = pi_real(m)
if mod(m, 2) == 0
    Qm = 1/sqrt(2) * [
        eye(m/2),         1j * eye(m/2)
        flip(eye(m/2)),   -1j * flip(eye(m/2))
        ];
else
    Qm = 1/sqrt(2) * [
        eye((m-1)/2),       zeros((m-1)/2, 1), 1j * eye((m-1)/2)
        zeros(1, (m-1)/2),  sqrt(2),           zeros(1, (m-1)/2)
        flip(eye((m-1)/2)), zeros((m-1)/2, 1), -1j * flip(eye((m-1)/2))
        ];
end

end