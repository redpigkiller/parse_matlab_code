function [r_eZF_Mc, r_eMMSE_Mc, r_eCORR_Mc] = exp_one_Mc_ber(i_Mc, sim_Mc_info, DEBUG)
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
r_ber_Mc = zeros(n_Q, n_snr);
r_eZF_Mc = cell(n_Q, n_snr);
r_eMMSE_Mc = cell(n_Q, n_snr);
r_eCORR_Mc = cell(n_Q, n_snr);

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
		r_ber_trial = zeros(n_trial_q, n_data_block);
		r_ZF_trial = cell(n_trial_q, n_data_block);
		r_MMSE_trial = cell(n_trial_q, n_data_block);
		r_CORR_trial = cell(n_trial_q, n_data_block);

		for i_trial_q = 1 : n_trial_q
			%% unpack eps_q_trial & h_q_trial
			eps_q_trial = eps_q(i_trial_q, :).';
			h_q_trial = permute(h_q(i_trial_q, :, :), [2, 3, 1]);
% size(eps_q)
% size(eps_hat_q)
% i_trial_q
			eps_hat_q_trial = eps_hat_q(i_trial_q, :).';
			h_hat_q_trial = permute(h_hat_q(i_trial_q, :, :), [2, 3, 1]);


chl_mtx_est = gen_channel_matrix(eps_hat_q_trial, h_hat_q_trial, N);
chl_mtx_cor = gen_channel_matrix(eps_q_trial, h_q_trial, N);

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

% choose equalizer
	
%             %%%%%%%%%%%%%%%%%%%%%%%
% MMSE
%             %%%%%%%%%%%%%%%%%%%%%%%
N0_yp = 0;
for i_block = 1 : n_data_block
    N0_yp = N0_yp + abs(yp_cell{i_block}' * yp_cell{i_block} - norm(chl_mtx_est, 'fro')^2) / N;
end
N0_yp = N0_yp / n_data_block;

MMSE_mtx_est = chl_mtx_est' * inverse(full(chl_mtx_est * chl_mtx_est' + N0_yp * speye(N)));
MMSE_mtx_est = ifft(fft(MMSE_mtx_est), N, 2);

MMSE_mtx_cor = chl_mtx_cor' * inverse(full(chl_mtx_cor * chl_mtx_cor' + 10^(-snr(i_snr)/10) * speye(N)));
MMSE_mtx_cor = ifft(fft(MMSE_mtx_cor), N, 2);
				

chl_mtx_est = ifft(fft(full(chl_mtx_est)), N, 2);
chl_mtx_cor = ifft(fft(full(chl_mtx_cor)), N, 2);

%             %%%%%%%%%%%%%%%%%%%%%%%
% ZF
%             %%%%%%%%%%%%%%%%%%%%%%%
ZF_mtx_est = inverse(chl_mtx_est);
ZF_mtx_cor = inverse(chl_mtx_cor);

            %% recover Tx data
			for i_block = 1 : n_data_block
                %% unpack s_block & yp
				s_block = s_data(i_block, :).';
                yp = yp_cell{i_block};

% 2-2. estimate Tx data (frequency domain)
s_ZF_raw_est = ZF_mtx_est * yp;
s_ZF_raw_idl = ZF_mtx_cor * yp;
s_MMSE_raw_est = MMSE_mtx_est * yp;
s_MMSE_raw_idl = MMSE_mtx_cor * yp;

% 3. demodulation (hard decision)
ss_ZF_est = nrSymbolDemodulate(s_ZF_raw_est, 'QPSK', 'DecisionType', 'Hard');
ss_ZF_idl = nrSymbolDemodulate(s_ZF_raw_idl, 'QPSK', 'DecisionType', 'Hard');
ss_MMSE_est = nrSymbolDemodulate(s_MMSE_raw_est, 'QPSK', 'DecisionType', 'Hard');
ss_MMSE_idl = nrSymbolDemodulate(s_MMSE_raw_idl, 'QPSK', 'DecisionType', 'Hard');

% estimate symbol
s_ZF_est = nrSymbolModulate(ss_ZF_est, 'QPSK');
s_ZF_idl = nrSymbolModulate(ss_ZF_idl, 'QPSK');
s_MMSE_est = nrSymbolModulate(ss_MMSE_est, 'QPSK');
s_MMSE_idl = nrSymbolModulate(ss_MMSE_idl, 'QPSK');
% correct symbol
s_cor = nrSymbolModulate(s_block, 'QPSK');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0
e0_ZF_est = norm(s_ZF_raw_est - s_ZF_est).^2;
e0_ZF_idl = norm(s_ZF_raw_idl - s_ZF_idl).^2;
e0_MMSE_est = norm(s_MMSE_raw_est - s_MMSE_est).^2;
e0_MMSE_idl = norm(s_MMSE_raw_idl - s_MMSE_idl).^2;

% 1
e1_ZF_est = yp - chl_mtx_est * s_ZF_est;
e1_MMSE_est = yp - chl_mtx_est * s_MMSE_est;

% 2

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


				ber = biterr(s_block, ss_MMSE_est) / numel(s_block);
% check
aa = find(s_block~=ss_MMSE_est);
bb = biterr(s_block, ss_MMSE_est)  ;   
if length(aa)~= bb
    error("11");
end
aa2 = unique(ceil(aa./2));
bb2 = find(s_MMSE_est~=s_cor);
if length(aa2)~= length(bb2)
    error("11");
end
aa3 = 1 : N;
aa3(aa2) = [];
bb3 = find(s_MMSE_est==s_cor);
if length(aa3)~= length(bb3)
    error("11");
end

aa = find(s_block~=ss_ZF_est);
bb = biterr(s_block, ss_ZF_est)  ;      
if length(aa)~= bb
    error("11");
end
aa2 = unique(ceil(aa./2));
bb2 = find(s_ZF_est~=s_cor);
if length(aa2)~= length(bb2)
    error("11");
end
aa3 = 1 : N;
aa3(aa2) = [];
bb3 = find(s_ZF_est==s_cor);
if length(aa3)~= length(bb3)
    error("11");
end


% ZF
err_idx = find(s_ZF_est~=s_cor);
cor_idx = find(s_ZF_est==s_cor);
if length(err_idx) + length(cor_idx) ~= N
    error("11");
end

eVal = abs(chl_mtx_est' * e1_ZF_est).^2;
errVal_ZF = eVal(err_idx);
corVal_ZF = eVal(cor_idx);

eVal = abs(s_ZF_est - s_ZF_raw_est).^2;
errVal2_ZF = eVal(err_idx);
corVal2_ZF = eVal(cor_idx);

% MMSE
err_idx = find(s_MMSE_est~=s_cor);
cor_idx = find(s_MMSE_est==s_cor);
if length(err_idx) + length(cor_idx) ~= N
    error("11");
end

eVal = abs(chl_mtx_est' * e1_MMSE_est).^2;
errVal_MMSE = eVal(err_idx);
corVal_MMSE = eVal(cor_idx);

eVal = abs(s_MMSE_est - s_MMSE_raw_est).^2;
errVal2_MMSE = eVal(err_idx);
corVal2_MMSE = eVal(cor_idx);

				if ber > 0.5
					ber = 0.5;
				end

				% save decoded data
				r_ber_trial(i_trial_q, i_block) = ber;

% ZF
errVal = struct;

errVal.e0_est = e0_ZF_est;
errVal.e0_idl = e0_ZF_idl;

errVal.e1_err_pts = errVal_ZF;
errVal.e1_cor_pts = corVal_ZF;

errVal.e2_err_pts = errVal2_ZF;
errVal.e2_cor_pts = corVal2_ZF;

r_ZF_trial{i_trial_q, i_block} = errVal;

% MMSE
errVal = struct;

errVal.e0_est = e0_MMSE_est;
errVal.e0_idl = e0_MMSE_idl;

errVal.e1_err_pts = errVal_MMSE;
errVal.e1_cor_pts = corVal_MMSE;

errVal.e2_err_pts = errVal2_MMSE;
errVal.e2_cor_pts = corVal2_MMSE;

r_MMSE_trial{i_trial_q, i_block} = errVal;

% r_CORR_trial{i_trial_q, i_block} = norm(e_corr)^2;

			end
		end

		%% save record
		r_ber_Mc(i_Q, i_snr) = sum(r_ber_trial, 'all') / numel(r_ber_trial);

r_eZF_Mc{i_Q, i_snr} = r_ZF_trial;
r_eMMSE_Mc{i_Q, i_snr} = r_MMSE_trial;
r_eCORR_Mc{i_Q, i_snr} = r_CORR_trial;
	end
end
end
