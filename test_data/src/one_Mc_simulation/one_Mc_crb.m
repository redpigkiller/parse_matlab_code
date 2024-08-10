function [r_dfo_crb, r_aoa_crb] = one_Mc_crb(i_Mc, sim_Mc_info, skip_Q, skip_snr)
if nargin < 3
    skip_Q      = [];
    skip_snr    = [];
elseif nargin < 4
    skip_snr    = [];
end

% generate more snapshots
snapshots_factor = 1;

%% unpack variables
snr         = sim_Mc_info.snr;
Q           = sim_Mc_info.Q;
Nr          = sim_Mc_info.Nr;
P           = sim_Mc_info.P;
N           = sim_Mc_info.config.N;
SP          = sim_Mc_info.config.SP;

% path of data
dir_gendata = sim_Mc_info.dir_gendata;

% index
snr_idx = sim_Mc_info.snr_idx;
Q_idx   = sim_Mc_info.Q_idx;

%% auxiliary variables
M       = N / P;
n_snr   = length(snr);
n_Q     = length(Q);

eps_q_max    = sim_Mc_info.config.fc / sim_Mc_info.config.fs * ...
    (sim_Mc_info.config.v / 3.6) / physconst('LightSpeed');

%% load data
chl_data = load(sprintf("%s_%03d.mat", dir_gendata, i_Mc));

%% unpack data
dfoVal = chl_data.dfoVal;
aoaVal = chl_data.aoaVal;
chlVal = chl_data.chlVal;

%% recording variables
r_dfo_crb = zeros(n_Q, n_snr);
r_aoa_crb = zeros(n_Q, n_snr);

for i_snr = 1 : n_snr
	%% pre-generate signal and noise (trainging data)
	% training data
    x_training = cell(snapshots_factor, 1);
    for ii = 1 : snapshots_factor
        s = randi([0 1], 2*M, 1);
	    x = nrSymbolModulate(s, 'QPSK');
	    x = sqrt(M) * ifft(x);

        x_training{ii} = repmat(x, P, 1);
    end

	for i_Q = 1 : n_Q
        if ismember(snr(i_snr), skip_snr) && ismember(Q(i_Q), skip_Q)
            continue;
        end

		%% unpack data
		dfo_q = dfoVal{Q_idx(i_Q), snr_idx(i_snr)};
		aoa_q = aoaVal{Q_idx(i_Q), snr_idx(i_snr)};
		chl_q = chlVal{Q_idx(i_Q), snr_idx(i_snr)};

		%% declare record var. for trial
		r_dfo_crb_trial = 0;
		r_aoa_crb_trial = 0;

        n_trial_q = size(dfo_q, 1);
        if (n_trial_q ~= size(aoa_q, 1) || n_trial_q ~= size(chl_q, 1))
            error("unmatched # of trials");
        end

        for i_trial_q = 1 : n_trial_q
			%% unpack eps_q & theta_q & h_q
			eps_q   = dfo_q(i_trial_q, :).';
			theta_q = aoa_q(i_trial_q, :).';
			h_q     = permute(chl_q(i_trial_q, :, :), [2, 3, 1]);

            if Nr == 1
                theta_q = zeros(Q(i_Q), 1);
            end
            % 1. compute covariance matrix
%             Rx = zeros(Q(i_Q));
%             for ii = 1 : snapshots_factor
%                 X = get_src_mtx(x_training{ii}, eps_q, h_q, Q(i_Q), M, N);
%                 for jj = 1 : M
%                     Rx  = Rx + X(:, jj) * X(:, jj)';
%                 end
%             end
%             Rx = 1/(snapshots_factor*M) * Rx;

            % theoretic Rx
            B = cell(Q(i_Q), 1);
            for jj = 1 : Q(i_Q)
                B{jj} = exp(2j*pi*eps_q(jj)) * circulant(h_q(jj, :).', M);
            end
            Rx = zeros(Q(i_Q));
            for jj = 1 : Q(i_Q)
                for kk = jj : Q(i_Q)
                    Rx(jj, kk) = 1/M * trace(B{jj} * B{kk}');
                end
                for kk = 1 : jj - 1
                    Rx(jj, kk) = conj(Rx(kk, jj));
                end
            end

            % 2. noise variance
            noise_var   = 10^(-snr(i_snr)/10);

            % 3. array steering matrix and its derivatives
            if Nr == 1
                A   = exp(1j * 2 * pi * eps_q * M * (0:P-1)).';
                D   = (1j*2*pi*M*(0:P-1).') .* A;

            else
                a_q1    = exp(-1j * 2 * pi * SP * theta_q * (0:Nr-1)).';
                a_q2    = exp(1j * 2 * pi * eps_q * M * (0:P-1)).';
                A       = kr(a_q1, a_q2);

                D1  = kr((-1j*2*pi*SP*(0:Nr-1)).' .* a_q1, a_q2);
                D2  = kr(a_q1, (1j*2*pi*M*(0:P-1)).' .* a_q2);
                D   = [D1 D2];
            end

            % 4. compute CRB matrix
            if Nr == 1
                CRB_mtx = get_crb_1d_UMA(A, D, Rx, noise_var, M);
            else
                CRB_mtx = get_crb_2d_UMA(A, D, Rx, noise_var, M);
            end

            % 5. average CRB
            crb = diag(CRB_mtx);
            % upper bound of CRB of DFO & AoA
            bnd_dfo = (0.5/M)^2/3;
            % bnd_dfo = eps_q_max^2/2;
            bnd_aoa = 1/3;

            % obtain crb
            if Nr == 1
                crb_dfo = crb;
            else
                crb_aoa = crb(1 : Q(i_Q));
                crb_dfo = crb(Q(i_Q)+1 : end);
            end

            % bound
            crb_dfo(crb_dfo > bnd_dfo)  = bnd_dfo;
            crb_dfo(crb_dfo < -bnd_dfo) = -bnd_dfo;
            if Nr > 1
                crb_aoa(crb_aoa > bnd_aoa)  = bnd_aoa;
                crb_aoa(crb_aoa < -bnd_aoa) = -bnd_aoa;
            end

            % average
            if Nr == 1
                crb_aoa = 0;
                crb_dfo = mean(abs(crb_dfo));
            else
                crb_aoa = mean(abs(crb_aoa));
                crb_dfo = mean(abs(crb_dfo));
            end

            %% record CRB
			r_dfo_crb_trial = r_dfo_crb_trial + crb_dfo;
			r_aoa_crb_trial = r_aoa_crb_trial + crb_aoa;

        end

		%% save record
		r_dfo_crb(i_Q, i_snr) = r_dfo_crb_trial / n_trial_q;
		r_aoa_crb(i_Q, i_snr) = r_aoa_crb_trial / n_trial_q;
	end
end

end

%% functions
function CRB_mtx = get_crb_2d_UMA(A, D, P, noise_var, snapshots)
% variables name mapping:
%   1. Q        <--->   n
%   2. Nr * P   <--->   m
%
n   = size(P, 1);
m   = size(A, 1);

R = A * P * A' + noise_var * eye(m);
H = D' * (eye(m) - A / (A' * A) * A') * D;
U = P * A' / R * A * P;
CRB_mtx = real(H .* (kron(ones(2), U).'));
CRB_mtx = noise_var / 2 / snapshots * ...
    (eye(2*n) / CRB_mtx);
end

function CRB_mtx = get_crb_1d_UMA(A, D, P, noise_var, snapshots)
% variables name mapping:
%   1. Q    <--->   n
%   2. P    <--->   m
%
n   = size(P, 1);
m   = size(A, 1);

R = A * P * A' + noise_var * eye(m);
H = D' * (eye(m) - A / (A' * A) * A') * D;
U = P * A' / R * A * P;
CRB_mtx = real(H .* (U.'));
CRB_mtx = noise_var / 2 / snapshots * ...
    (eye(n) / CRB_mtx);

end

function X = get_src_mtx(x, eps_q, h_q, Q, M, N)
X = zeros(Q, M);

for i_Q = 1 : Q
    h_q_vec = h_q(i_Q, :);
    yp = cconv(x.', h_q_vec, N);
    yp = yp(1 : M);
    X(i_Q, :) = yp .* exp(2j * pi * eps_q(i_Q) * (0:M-1));
end
end
