function chl_mtx = gen_channel_matrix(dfo, aoa, chl, N, Q, Nr, SP, algo_flag, is_sparse)
%GEN_CHANNEL_MATRIX Generate the channel matrix.
%
%   See also CIRCULANT.
%

if algo_flag.assume_no_dfo == 1
    dfo = zeros(Q, 1);
end
if algo_flag.assume_no_aoa == 1
    aoa = zeros(Q, 1);
end

if nargin < 9
    is_sparse = 1;
end

%% check
assert(Q == length(dfo), "unmatch number of paths");
if Nr > 1
    assert(Q == length(aoa), "unmatched Q");
end
assert(Q == size(chl, 1), "unmatched Q");

%% generate channel matrix
% generate array manifold matrix
if Nr == 1
    A_hat = exp(-1j * 2 * pi * SP * zeros(1, Q, 'like', dfo));
else
    A_hat = exp(-1j * 2 * pi * SP * aoa * (0:Nr-1)).';
end

% generate DFO and channel matrix
if is_sparse == 0
    B_hat = zeros(N*Q, N, 'like', dfo);

    chl = [chl zeros(Q, N - size(chl, 2), 'like', dfo)];
    
    for i_Q = 1 : Q
	    DFO_q   = diag(exp(1j * 2 * pi * dfo(i_Q) * (0 : N - 1)));
	    H_q     = circulant(chl(i_Q, :).', N);
	    B_hat((i_Q-1)*N+1 : i_Q*N, :) = DFO_q * H_q;
    end

else
    B_hat = cell(Q, 1);

    for i_Q = 1 : Q
	    DFO_q = spdiags(exp(1j * 2 * pi * dfo(i_Q) * (0 : N - 1)).', 0, N, N);

        % construct sparse(chl(i_Q, :)) with size N x 1
        [chl_i, chl_j, chl_v] = find(chl(i_Q, :).');
        sph = sparse(chl_i, chl_j, chl_v, N, 1);

	    H_q = circulant(sph, N);

        % store the indices and values
        B_hat{i_Q} = DFO_q * H_q;
        
    end

    B_hat = vertcat(B_hat{:});

end

% generate completed matrix
chl_mtx = kron(A_hat, speye(N)) * B_hat;

end