function [r_MSE_dfo, r_MSE_aoa, r_MSE_chl, s_est_dfo, s_est_aoa, s_est_chl] = one_Mc_est(i_Mc, sim_Mc_info, skip_Q, skip_snr)
%ONE_MC_EST The function simulates for one Monte-Carlo with specfied i_Mc
%   sim_Mc_info: the data & information needed to perform simulation
%

if nargin < 3
    skip_Q      = [];
    skip_snr    = [];
elseif nargin < 4
    skip_snr    = [];
end

%% unpack variables
snr         = sim_Mc_info.snr;
Q           = sim_Mc_info.Q;
Nr          = sim_Mc_info.Nr;
P           = sim_Mc_info.P;
N           = sim_Mc_info.config.N;
SP          = sim_Mc_info.config.SP;
Lch         = sim_Mc_info.config.Lch;
add_CP_len  = iif(sim_Mc_info.glob_flag.enable_CP_sim, sim_Mc_info.config.Lcp, 0);

% flags
flag_use_gpu    = sim_Mc_info.glob_flag.use_gpu;
flag_add_noise  = sim_Mc_info.glob_flag.add_noise;
algo_flag       = sim_Mc_info.algo_flag;

% path of data
dir_gendata = sim_Mc_info.dir_gendata;

% index
snr_idx = sim_Mc_info.snr_idx;
Q_idx   = sim_Mc_info.Q_idx;

%% auxiliary variables
M                   = N / P;
n_snr               = length(snr);
n_Q                 = length(Q);

if flag_use_gpu
	p_type = zeros('gpuArray');
else
	p_type = zeros();
end

%% load data
chl_data = load(sprintf("%s_%03d.mat", dir_gendata, i_Mc));

%% unpack data
dfoVal = chl_data.dfoVal;
aoaVal = chl_data.aoaVal;
chlVal = chl_data.chlVal;

%% recording variables
r_MSE_dfo = zeros(n_Q, n_snr);
r_MSE_aoa = zeros(n_Q, n_snr);
r_MSE_chl = zeros(n_Q, n_snr);
s_est_dfo = cell(n_Q, n_snr);
s_est_aoa = cell(n_Q, n_snr);
s_est_chl = cell(n_Q, n_snr);

for i_snr = 1 : n_snr
	%% pre-generate signal and noise (trainging data)
	% 1. training data
	s_training = randi([0 1], 2*M, 1);
	x_training = nrSymbolModulate(s_training, 'QPSK');
	x_training = sqrt(M) * ifft(x_training);
	if flag_use_gpu;	x_training = gpuArray(x_training);	end
	x_training = repmat(x_training, P, 1);

	% 2. noise term
	sigma_w     = sqrt(10.^(-snr(i_snr)/ 10) / 2);
	n_training  =        sigma_w .* randn(N * Nr, 1, 'like', p_type) + ...
		            1j * sigma_w .* randn(N * Nr, 1, 'like', p_type);

	for i_Q = 1 : n_Q
        if ismember(snr(i_snr), skip_snr) && ismember(Q(i_Q), skip_Q)
            continue;
        end

		%% unpack data
		dfo_q = dfoVal{Q_idx(i_Q), snr_idx(i_snr)};
		aoa_q = aoaVal{Q_idx(i_Q), snr_idx(i_snr)};
		chl_q = chlVal{Q_idx(i_Q), snr_idx(i_snr)};

		if flag_use_gpu;	dfo_q = gpuArray(dfo_q);	end
		if flag_use_gpu;	aoa_q = gpuArray(aoa_q);	end
		if flag_use_gpu;	chl_q = gpuArray(chl_q);	end

		%% declare record var. for trial
		r_MSE_dfo_trial = 0;
		r_MSE_aoa_trial = 0;
		r_MSE_chl_trial = 0;

        n_trial_q = size(dfo_q, 1);
        if (n_trial_q ~= size(aoa_q, 1) || n_trial_q ~= size(chl_q, 1))
            error("unmatched # of trials");
        end

		s_est_dfo_trial = zeros(n_trial_q, Q(i_Q));
		s_est_aoa_trial = zeros(n_trial_q, Q(i_Q));
		s_est_chl_trial = zeros(n_trial_q, Q(i_Q), Lch+1);

		for i_trial_q = 1 : n_trial_q
			%% unpack eps_q & theta_q & h_q
			eps_q   = dfo_q(i_trial_q, :).';
			theta_q = aoa_q(i_trial_q, :).';
			h_q     = permute(chl_q(i_trial_q, :, :), [2, 3, 1]);

			% 1. compute Rx data (frequency domain)
			rp = dfo_channel_ula(...
				x_training, n_training, ...
				eps_q, theta_q, h_q, ...
				N, Q(i_Q), Nr, SP, ...
                flag_add_noise, add_CP_len ...
			);

			% 2. estimate channel
            [eps_hat_q, theta_hat_q, h_hat_q] = estimate_channel(...
				rp, ...
				x_training, ...
				Q(i_Q), P, Nr, ...
				sim_Mc_info.config, ...
				algo_flag ...
			);

			%% evaluation MSE DFO and MSE channel
            % 3. pair the correction DFO & AOA and estimated DFO
			[eps_hat_q, pair_Ind]   = arrpair(eps_q, eps_hat_q, 0.5);
            theta_hat_q             = theta_hat_q(pair_Ind);
            h_hat_q                 = h_hat_q(pair_Ind, :);

			% 4. save estimated DFO and channel
			s_est_dfo_trial(i_trial_q, :)       = eps_hat_q;
			s_est_aoa_trial(i_trial_q, :)       = theta_hat_q;
			s_est_chl_trial(i_trial_q, :, :)    = h_hat_q(:, 1 : Lch + 1);

            % 5. calculate MSE
			MSE_eps         = norm(circdist(eps_hat_q, eps_q, 0.5, true))^2 / Q(i_Q);
			MSE_theta       = norm(circdist(theta_hat_q, theta_q, 1, true))^2 / Q(i_Q);
			mean_MSE_h_q    = norm(h_hat_q(:, 1 : Lch + 1) ...
				                 - h_q(:, 1 : Lch + 1), 'fro')^2 / Q(i_Q);

			%% record MSE
			r_MSE_dfo_trial = r_MSE_dfo_trial + MSE_eps;
			r_MSE_aoa_trial = r_MSE_aoa_trial + MSE_theta;
			r_MSE_chl_trial = r_MSE_chl_trial + mean_MSE_h_q;

		end
		%% save estimated DFO & AOA and channel
		s_est_dfo{i_Q, i_snr} = s_est_dfo_trial;
		s_est_aoa{i_Q, i_snr} = s_est_aoa_trial;
		s_est_chl{i_Q, i_snr} = s_est_chl_trial;

		%% save record
		r_MSE_dfo(i_Q, i_snr) = r_MSE_dfo_trial / n_trial_q;
		r_MSE_aoa(i_Q, i_snr) = r_MSE_aoa_trial / n_trial_q;
		r_MSE_chl(i_Q, i_snr) = r_MSE_chl_trial / n_trial_q;
	end
end

end