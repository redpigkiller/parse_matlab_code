function r_ber_Mc = one_Mc_ber(i_Mc, sim_Mc_info)
%ONE_MC_BER The function simulates for one Monte-Carlo with specfied i_Mc
%   sim_Mc_info: the data & information needed to perform simulation

% When Nr = 1, using ZF receiver, it might incur the rank deficient
% problem. Since the ZF receiver would calculate the LS solution, the
% solution would be ss_hat = pinv(chl_mtx) * yp.
%
% Note: for m x n matrix A and a linear system A * x = b
% + For overdetermined system (m > n), there might be no solution:
%   - mldivide, pinv, and lsqminnorm find the LS solution, i.e. minimize
%     norm(A*x - b)
%   - 
% + For underdetermined system (m < n), there are many solution:
%   - mldivide: find the solution that minimize norm(A*x - b) and 0-norm of
%     x (L0 norm of x), i.e. it also maximizes the number of zeros of x
%   - pinv: find the LS solution that minimize norm(A*x - b)
%   - lsqminnorm: find the solution that minimize norm(A*x - b) and 2-norm 
%     of x (L2 norm of x)
%

%% unpack variables
snr             = sim_Mc_info.snr;
Q               = sim_Mc_info.Q;
Nr              = sim_Mc_info.Nr;
N               = sim_Mc_info.config.N;
SP              = sim_Mc_info.config.SP;
add_CP_len      = iif(sim_Mc_info.glob_flag.enable_CP_sim, sim_Mc_info.config.Lcp, 0);
n_data_block    = sim_Mc_info.config.data_block;

% flags
flag_use_gpu    = sim_Mc_info.glob_flag.use_gpu;
flag_add_noise  = sim_Mc_info.glob_flag.add_noise;
algo_flag       = sim_Mc_info.algo_flag;

% path of data
dir_gendata = sim_Mc_info.dir_gendata;
dir_estdata = sim_Mc_info.dir_estdata;

% index
snr_gendata_idx = sim_Mc_info.snr_gendata_idx;
snr_estchal_idx = sim_Mc_info.snr_estchal_idx;
Q_gendata_idx   = sim_Mc_info.Q_gendata_idx;
Q_estchal_idx   = sim_Mc_info.Q_estchal_idx;

%% auxiliary variables
n_snr   = length(snr);
n_Q     = length(Q);

if flag_use_gpu
	p_type = zeros('gpuArray');
else
	p_type = zeros();
end

%% load data
chl_data    = load(sprintf("%s_%03d.mat", dir_gendata, i_Mc));
est_data    = load(sprintf("%s_est_%03d.est.mat", dir_estdata, i_Mc));

%% unpack data
dfoVal = chl_data.dfoVal;
aoaVal = chl_data.aoaVal;
chlVal = chl_data.chlVal;

est_dfoVal = est_data.est_dfoVal;
est_aoaVal = est_data.est_aoaVal;
est_chlVal = est_data.est_chlVal;

%% recording variables
r_ber_Mc = zeros(n_Q, n_snr);

for i_snr = 1 : n_snr
	%% pre-generate signal and noise (trainging data)
	% 1. transmitted data
	s_data = randi([0, 1], n_data_block, 2 * N);

	% 2. noise term
	sigma_w = sqrt(10.^(-snr(i_snr)/ 10) / 2);
	n_data  =        sigma_w .* randn(n_data_block, N * Nr, 'like', p_type) + ...
	            1j * sigma_w .* randn(n_data_block, N * Nr, 'like', p_type);

	for i_Q = 1 : n_Q
		%% unpack data
		dfo_q   = dfoVal{Q_gendata_idx(i_Q), snr_gendata_idx(i_snr)};
		aoa_q   = aoaVal{Q_gendata_idx(i_Q), snr_gendata_idx(i_snr)};
		chl_q   = chlVal{Q_gendata_idx(i_Q), snr_gendata_idx(i_snr)};

		dfo_hat_q = est_dfoVal{Q_estchal_idx(i_Q), snr_estchal_idx(i_snr)};
		aoa_hat_q = est_aoaVal{Q_estchal_idx(i_Q), snr_estchal_idx(i_snr)};
		chl_hat_q = est_chlVal{Q_estchal_idx(i_Q), snr_estchal_idx(i_snr)};

		if flag_use_gpu;	dfo_q = gpuArray(dfo_q);	end
		if flag_use_gpu;	aoa_q = gpuArray(aoa_q);	end
		if flag_use_gpu;	chl_q = gpuArray(chl_q);	end

		if flag_use_gpu;	dfo_hat_q = gpuArray(dfo_hat_q);	end
		if flag_use_gpu;	aoa_hat_q = gpuArray(aoa_hat_q);	end
		if flag_use_gpu;	chl_hat_q = gpuArray(chl_hat_q);	end

		%% declare record var. for trial

        n_trial_q = size(dfo_q, 1);
        if (n_trial_q ~= size(aoa_q, 1) || ...
            n_trial_q ~= size(chl_q, 1) || ...
            n_trial_q ~= size(dfo_hat_q, 1) || ...
            n_trial_q ~= size(aoa_hat_q, 1) || ...
            n_trial_q ~= size(chl_hat_q, 1) ...
        )
            error("unmatched # of trials");
        end

		r_ber_trial = zeros(n_trial_q, n_data_block);

		for i_trial_q = 1 : n_trial_q
			%% unpack
			eps_q   = dfo_q(i_trial_q, :).';
			theta_q = aoa_q(i_trial_q, :).';
			h_q     = permute(chl_q(i_trial_q, :, :), [2, 3, 1]);

			eps_hat_q   = dfo_hat_q(i_trial_q, :).';
			theta_hat_q = aoa_hat_q(i_trial_q, :).';
			h_hat_q     = permute(chl_hat_q(i_trial_q, :, :), [2, 3, 1]);

            switch (sim_Mc_info.glob_flag.gen_chl_mtx_flags)
                case 0
                    chl_mtx = gen_channel_matrix(eps_hat_q, theta_hat_q, h_hat_q, N, Q(i_Q), Nr, SP, algo_flag);
                case 1
                    chl_mtx = gen_channel_matrix(eps_q, theta_hat_q, h_hat_q, N, Q(i_Q), Nr, SP, algo_flag);
                case 2
                    chl_mtx = gen_channel_matrix(eps_hat_q, theta_q, h_hat_q, N, Q(i_Q), Nr, SP, algo_flag);
                case 3
                    chl_mtx = gen_channel_matrix(eps_q, theta_q, h_hat_q, N, Q(i_Q), Nr, SP, algo_flag);
                case 4
                    chl_mtx = gen_channel_matrix(eps_hat_q, theta_hat_q, h_q, N, Q(i_Q), Nr, SP, algo_flag);
                case 5
                    chl_mtx = gen_channel_matrix(eps_q, theta_hat_q, h_q, N, Q(i_Q), Nr, SP, algo_flag);
                case 6
                    chl_mtx = gen_channel_matrix(eps_hat_q, theta_q, h_q, N, Q(i_Q), Nr, SP, algo_flag);
                case 7
                    chl_mtx = gen_channel_matrix(eps_q, theta_q, h_q, N, Q(i_Q), Nr, SP, algo_flag);
                otherwise
                    error('error gen_chl_mtx_flags');
            end

            % pre-compute all yp
            yp_cell = cell(n_data_block, 1);
            for i_block = 1 : n_data_block
                %% unpack s_block & n_block
				s_block = s_data(i_block, :).';
				n_block = n_data(i_block, :).';

				%% generate Tx data
				x_transmitted = nrSymbolModulate(s_block, 'QPSK');
				if flag_use_gpu;	x_transmitted = gpuArray(x_transmitted);	end
				x_transmitted = sqrt(N) * ifft(x_transmitted);

				n_transmitted = n_block;
				if flag_use_gpu;	n_transmitted = gpuArray(n_transmitted);	end

				% 1. compute Rx data (frequency domain)
                rp = dfo_channel_ula(...
				    x_transmitted, n_transmitted, ...
				    eps_q, theta_q, h_q, ...
				    N, Q(i_Q), Nr, SP, ...
                    flag_add_noise, add_CP_len ...
			    );

                yp_cell{i_block} = rp;
            end

            % choose equalizer:
            flag_illcond = false;
            if strcmpi(algo_flag.equalizer, 'ZF')
                pchl_mtx = decomposition(chl_mtx);

                if isIllConditioned(pchl_mtx)
                    flag_illcond = true;
                    pchl_mtx = pinv(full(chl_mtx));
                end

			elseif strcmpi(algo_flag.equalizer, 'MMSE')
				if algo_flag.equ_opt.use_correct_snr == 1
                    pchl_mtx = chl_mtx' / (full(chl_mtx * chl_mtx' + 10^(-snr(i_snr)/10) * speye(N)));
                    
                elseif algo_flag.equ_opt.compute_snr_each_datablock == 0
                    N0_yp = 0;
                    for i_block = 1 : n_data_block
                        N0_yp = N0_yp + abs(yp_cell{i_block}' * yp_cell{i_block} - norm(chl_mtx, 'fro')^2) / N / Nr;
                    end
                    N0_yp = N0_yp / n_data_block;
                    
                    pchl_mtx = chl_mtx' / (full(chl_mtx * chl_mtx' + N0_yp * speye(N*Nr)));
				end
            end

            %% recover Tx data
			for i_block = 1 : n_data_block
                %% unpack s_block & yp
				s_block = s_data(i_block, :).';
                yp      = yp_cell{i_block};

				% 2-1. if use MMSE equalizer, need to estimate noise power
                if strcmpi(algo_flag.equalizer, 'MMSE') && algo_flag.equ_opt.compute_snr_each_datablock == 1
					N0_yp = abs(yp' * yp - norm(chl_mtx, 'fro')^2) / N;
					pchl_mtx = chl_mtx' / (full(chl_mtx * chl_mtx' + N0_yp * speye(N*Nr)));
                end

				% 2-2. estimate Tx data (frequency domain)
                if strcmpi(algo_flag.equalizer, 'ZF')
                    if ~flag_illcond
                        ss_hat = pchl_mtx \ yp;
                    else
                        ss_hat = pchl_mtx * yp;
                    end

			    elseif strcmpi(algo_flag.equalizer, 'MMSE')
				    ss_hat = pchl_mtx * yp;
                end

                ss_hat = 1/sqrt(N) * fft(ss_hat);

                % 3. post-processing (demodulation)
                if strcmpi(algo_flag.post, 'NONE')
				    %  (hard decision)
			        s_hat = nrSymbolDemodulate(ss_hat, 'QPSK', 'DecisionType', 'Hard');

                elseif strcmpi(algo_flag.post, 'ML_SEARCH')
                    s_hat = ml_search(yp, ss_hat, chl_mtx_post, algo_flag.post_opt);

                end

				%% bit error rate
				ber = biterr(s_block, s_hat) / numel(s_block);

				if ber > 0.5
					ber = 0.5;
				end

				% save decoded data
				r_ber_trial(i_trial_q, i_block) = ber;
			end
		end

		%% save record
		r_ber_Mc(i_Q, i_snr) = sum(r_ber_trial, 'all') / numel(r_ber_trial);
	end
end
end
