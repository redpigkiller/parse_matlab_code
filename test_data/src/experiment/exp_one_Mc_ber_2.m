function [r_ber_Mc_n, r_ber_Mc_ml] = exp_one_Mc_ber_2(i_Mc, sim_Mc_info, DEBUG)
%SIM_BER_ONE_MC The function simulates for one Monte-Carlo with specfied
%i_Mc
%   sim_Mc_info: the data & information needed to perform simulation

if nargin < 3
    DEBUG = 0;
end

%% unpack variables
snr = sim_Mc_info.snr;
Q = sim_Mc_info.Q;
N = sim_Mc_info.config.N;
n_data_block = sim_Mc_info.config.data_block;

% flags
flag_use_gpu = sim_Mc_info.glob_flag.use_gpu;
flag_add_noise = sim_Mc_info.glob_flag.add_noise;
algo_flag = sim_Mc_info.algo_flag;

% path of data
dir_gendata = sim_Mc_info.dir_gendata;
dir_estdata = sim_Mc_info.dir_estdata;

% index
snr_gendata_idx = sim_Mc_info.snr_gendata_idx;
snr_estchal_idx = sim_Mc_info.snr_estchal_idx;
Q_gendata_idx = sim_Mc_info.Q_gendata_idx;
Q_estchal_idx = sim_Mc_info.Q_estchal_idx;

%% auxiliary variables
n_snr = length(snr);
n_Q = length(Q);

if flag_use_gpu
	p_type = zeros('gpuArray');
else
	p_type = zeros();
end

%% load data
[channel_data_idx, ii_Mc]= segment_load(i_Mc, 1);
channel_data = load(sprintf("%s_%03d.mat", dir_gendata, channel_data_idx));

[estimat_data_idx, jj_Mc]= segment_load(i_Mc, 1);
est_data = load(sprintf("%s_s_est_%03d.chl.mat", dir_estdata, estimat_data_idx));

%% unpack data
eps_Mc = channel_data.dfoVal;
if ndims(eps_Mc) >= 3 %#ok<ISMAT> 
	eps_Mc = eps_Mc(:, :, ii_Mc);
end
h_Mc = channel_data.chlVal;
if ndims(h_Mc) >= 3 %#ok<ISMAT> 
	h_Mc = h_Mc(:, :, ii_Mc);
end

eps_hat_Mc = est_data.est_dfoVal;
if ndims(eps_hat_Mc) >= 3 %#ok<ISMAT> 
	eps_hat_Mc = eps_hat_Mc(:, :, jj_Mc);
end
h_hat_Mc = est_data.est_chlVal;
if ndims(h_hat_Mc) >= 3 %#ok<ISMAT> 
	h_hat_Mc = h_hat_Mc(:, :, jj_Mc);
end

%% recording variables
r_ber_Mc_n = zeros(n_Q, n_snr);
r_ber_Mc_ml = zeros(n_Q, n_snr);

for i_snr = 1 : n_snr
	%% pre-generate signal and noise (trainging data)
	% 1. transmitted data
	s_data = randi([0, 1], n_data_block, 2 * N);

	% 2. noise term
	sigma_w = sqrt(10.^(-snr(i_snr)/ 10) / 2);
	n_data = zeros(n_data_block, N, 'like', p_type);
	for i_block = 1 : n_data_block
		n_data(i_block, :) = sigma_w .* randn(N, 1, 'like', p_type) + ...
			1j * sigma_w .* randn(N, 1, 'like', p_type);
	end

	for i_Q = 1 : n_Q
		%% unpack data
		eps_q = eps_Mc{Q_gendata_idx(i_Q), snr_gendata_idx(i_snr)};
		h_q = h_Mc{Q_gendata_idx(i_Q), snr_gendata_idx(i_snr)};

		eps_hat_q = eps_hat_Mc{Q_estchal_idx(i_Q), snr_estchal_idx(i_snr)};
		h_hat_q = h_hat_Mc{Q_estchal_idx(i_Q), snr_estchal_idx(i_snr)};

		if flag_use_gpu;	eps_q = gpuArray(eps_q);	end
		if flag_use_gpu;	h_q = gpuArray(h_q);	end

		if flag_use_gpu;	eps_hat_q = gpuArray(eps_hat_q);	end
		if flag_use_gpu;	h_hat_q = gpuArray(h_hat_q);	end

		%% declare record var. for trial
		n_trial_q = size(eps_q, 1);
		r_ber_trial_n = zeros(n_trial_q, n_data_block);
		r_ber_trial_ml = zeros(n_trial_q, n_data_block);

		for i_trial_q = 1 : n_trial_q
			%% unpack eps_q_trial & h_q_trial
			eps_q_trial = eps_q(i_trial_q, :).';
			h_q_trial = permute(h_q(i_trial_q, :, :), [2, 3, 1]);

			eps_hat_q_trial = eps_hat_q(i_trial_q, :).';
			h_hat_q_trial = permute(h_hat_q(i_trial_q, :, :), [2, 3, 1]);

			switch (sim_Mc_info.glob_flag.use_correct_dfo * 2 + sim_Mc_info.glob_flag.use_correct_chl)
				case 0
					chl_mtx = gen_channel_matrix(eps_hat_q_trial, h_hat_q_trial, N);
				case 1
					chl_mtx = gen_channel_matrix(eps_hat_q_trial, h_q_trial, N);
				case 2
					chl_mtx = gen_channel_matrix(eps_q_trial, h_hat_q_trial, N);
				case 3
					chl_mtx = gen_channel_matrix(eps_q_trial, h_q_trial, N);
			end

            if algo_flag.null_subcarrier > 0
				if algo_flag.position == 0
					chl_mtx = chl_mtx(:, algo_flag.null_subcarrier + 1 : ...
						end - algo_flag.null_subcarrier);
				else
					chl_mtx = chl_mtx(:, [1 : N/2 - algo_flag.null_subcarrier - 1, ...
						N/2 + algo_flag.null_subcarrier + 1 : N]);
				end
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

				if algo_flag.null_subcarrier > 0
					if algo_flag.position == 0
						x_transmitted(1 : algo_flag.null_subcarrier) = 0;
						x_transmitted(end - algo_flag.null_subcarrier + 1 : end) = 0;
					else
						x_transmitted(N/2 - algo_flag.null_subcarrier : N/2 + algo_flag.null_subcarrier) = 0;
					end
				end

				n_transmitted = n_block;
				if flag_use_gpu;	n_transmitted = gpuArray(n_transmitted);	end

				% 1. compute Rx data (frequency domain)
				rp = dfo_channel(...
					x_transmitted, ...
					n_transmitted, ...
					eps_q_trial, ...
					h_q_trial, ...
					flag_add_noise ...
				);

                yp_cell{i_block} = 1/sqrt(N) * fft(rp);
            end

			% choose equalizer:
            % for post-processing, store channel matrix in chl_mtx_post:
            %
            % if equalizer is ZF, chl_mtx is then not needed, discards it
            % after pchl_mtx generated.
            %
            % if equalizer is MMSE, the sparse chl_mtx is needed to
            % decrease computational cost, so keeps sparse chl_mtx and also
            % chl_mtx_post.

            if strcmpi(algo_flag.post, 'ML_SEARCH')
                chl_mtx_post = ifft(fft(full(chl_mtx)), N, 2);
            end
			if strcmpi(algo_flag.equalizer, 'ZF')
                if strcmpi(algo_flag.post, 'NONE')
                    chl_mtx = ifft(fft(full(chl_mtx)), N, 2);
                    pchl_mtx = inverse(chl_mtx);
                else
                    pchl_mtx = inverse(chl_mtx_post);
                end
				
			elseif strcmpi(algo_flag.equalizer, 'MMSE')
				if algo_flag.equ_opt.use_correct_snr == 1
                    pchl_mtx = chl_mtx' * inverse(full(chl_mtx * chl_mtx' + 10^(-snr(i_snr)/10) * speye(N)));
                    pchl_mtx = ifft(fft(pchl_mtx), N, 2);
                    
                elseif algo_flag.equ_opt.compute_snr_each_datablock == 0
                    N0_yp = 0;
                    for i_block = 1 : n_data_block
                        N0_yp = N0_yp + abs(yp_cell{i_block}' * yp_cell{i_block} - norm(chl_mtx, 'fro')^2) / N;
                    end
                    N0_yp = N0_yp / n_data_block;
                    
                    pchl_mtx = chl_mtx' * inverse(full(chl_mtx * chl_mtx' + N0_yp * speye(N)));
                    pchl_mtx = ifft(fft(pchl_mtx), N, 2);
				end
			end

            %% recover Tx data
			for i_block = 1 : n_data_block
                %% unpack s_block & yp
				s_block = s_data(i_block, :).';
                yp = yp_cell{i_block};

				% 2-1. if use MMSE equalizer, need to estimate noise power
				if strcmpi(algo_flag.equalizer, 'MMSE') && algo_flag.equ_opt.compute_snr_each_datablock == 1
					N0_yp = abs(yp' * yp - norm(chl_mtx, 'fro')^2) / N;
					pchl_mtx = chl_mtx' * inverse(full(chl_mtx * chl_mtx' + N0_yp * speye(N)));
                    pchl_mtx = ifft(fft(pchl_mtx), N, 2);
				end

				% 2-2. estimate Tx data (frequency domain)
				ss_hat = pchl_mtx * yp;

                % 3. post-processing (demodulation)
%                 if strcmpi(algo_flag.post, 'NONE')
				    %  (hard decision)
    s_hat_n = nrSymbolDemodulate(ss_hat, 'QPSK', 'DecisionType', 'Hard');

%                 elseif strcmpi(algo_flag.post, 'ML_SEARCH')
    s_hat_ml = ml_search(yp, ss_hat, chl_mtx_post, algo_flag.post_opt);

%                 end

				%% bit error rate
				if algo_flag.null_subcarrier > 0
					if algo_flag.position == 0
						s_block = s_block(2 * algo_flag.null_subcarrier + 1 : ...
							end - 2 * algo_flag.null_subcarrier);
					else
						s_block = s_block([1 : N - 2*algo_flag.null_subcarrier - 2, ...
							N + 2*algo_flag.null_subcarrier + 1 : 2*N]);
					end
                end


% find(s_hat_n~=s_block).'
% find(s_hat_ml~=s_block).'
				ber_n = biterr(s_block, s_hat_n) / numel(s_block);
				ber_ml = biterr(s_block, s_hat_ml) / numel(s_block);
%     if 2*N ~= numel(s_block)
%         gg
%     end
% ber_n
% ber_ml
% disp("2222222222222222222222222222222222222");
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s_hat1 = nrSymbolDemodulate(ss_hat, 'QPSK', 'DecisionType', 'Hard');
% ber1 = biterr(s_block, s_hat1) / numel(s_block);
% if (ber1 < ber)
%     xxxx1 = nrSymbolModulate(s_hat1, 'QPSK');
%     xxxx = nrSymbolModulate(s_hat, 'QPSK');
%     if (norm(yp-chl_mtx_post*xxxx1)^2 < norm(yp-chl_mtx_post*xxxx)^2)
%         error("111222222222");
%     end
% else
%     
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 				if ber > 0.5
% 					ber = 0.5;
% 				end

				% save decoded data
				r_ber_trial_n(i_trial_q, i_block) = ber_n;
				r_ber_trial_ml(i_trial_q, i_block) = ber_ml;
			end
		end

		%% save record
%         if n_trial_q * n_data_block ~= numel(r_ber_trial_n)
%             ggg
%         end
%         if n_trial_q * n_data_block ~= numel(r_ber_trial_ml)
%             gggg
%         end
        
		r_ber_Mc_n(i_Q, i_snr) = sum(r_ber_trial_n, 'all') / numel(r_ber_trial_n);
		r_ber_Mc_ml(i_Q, i_snr) = sum(r_ber_trial_ml, 'all') / numel(r_ber_trial_ml);
	end
end
end
