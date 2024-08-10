function [eps_hat_q, theta_hat_q, h_hat_q] = estimate_channel(yp, xp, Q, P, Nr, arg, flags)
%ESTIMATE_CHANNEL Estimate the channel based on yp, use the approch in flags
%   yp:     input signal matrix
%   xp:     (known) training sequence
%   Q:      # of paths
%   arg:    some configuration
%   flags:  flags of estimating channel
%

% auxiliary variables
M       = arg.N / P;
x       = xp(1 : M);
X_cir   = circulant(x, arg.Lch+1);

% partition yp into a (Nr x P) x M matrix
Y_hat = reshape(yp, M, []).';

%% %%%%%%%%%%%%%%%%%%%% estimate dfo %%%%%%%%%%%%%%%%%%%%
if flags.estimate_dfo == 1
	% estimate DFO using ESPRIT
    if Nr > 1
        if strcmpi(flags.use_dfo_aoa_algo, 'UNITARY_ESPRIT_2D')
            [eps_hat_q, theta_hat_q] = unitaryesprit_2d(Y_hat, P, Q, flags.dfo_aoa_opt);
            eps_hat_q   = eps_hat_q / M;
            theta_hat_q = -theta_hat_q / arg.SP;

            % separate eps_q and eps_hat_q (prevent bad inverse matrix)
            if flags.enable_sep ~= 0
                eps_q_max   = arg.fc / arg.fs * (arg.v / 3.6) / physconst('LightSpeed');
	            eps_hat_q   = arrsepoptim(eps_hat_q, arg.dfo_gap * eps_q_max);
	            theta_hat_q = arrsepoptim(theta_hat_q, arg.aoa_gap);
            end

	    else
		    error("wrong algorithm selected")
        end
        
    else

        if strcmpi(flags.use_dfo_aoa_algo, 'ESPRIT')
		    if flags.dfo_opt.minus_mean == 0
			    [eps_hat_q, ~] = esprit(Y_hat, Q, flags.dfo_opt);
		    else
			    [eps_hat_q, ~] = esprit(Y_hat - mean(Y_hat, 2), Q, ...
                    flags.dfo_opt);
		    end
		    eps_hat_q = eps_hat_q / M;
    
        elseif strcmpi(flags.use_dfo_aoa_algo, 'UNITARY_ESPRIT')
            [eps_hat_q, ~]  = unitaryesprit(Y_hat, Q, flags.dfo_opt);
            eps_hat_q       = eps_hat_q / M;
    
	    elseif strcmpi(flags.use_dfo_aoa_algo, 'ML_Estimation')
		    eps_hat_q = ml_estimation(Y_hat, Q, flags.dfo_opt, arg);
    
	    else
		    error("wrong algorithm selected")
        end

        % separate eps_q and eps_hat_q (prevent bad inverse matrix)
        if flags.enable_sep ~= 0
            eps_q_max   = arg.fc / arg.fs * (arg.v / 3.6) / physconst('LightSpeed');
            eps_hat_q   = arrsepoptim(eps_hat_q, arg.dfo_gap * eps_q_max); 
        end
        
        theta_hat_q = zeros(Q, 1, 'like', yp);
    end

else
	eps_hat_q   = zeros(Q, 1, 'like', yp);
	theta_hat_q = zeros(Q, 1, 'like', yp);
end

%% %%%%%%%%%%%%%%%%%%%% simulate channel %%%%%%%%%%%%%%%%%%%%
if flags.estimate_chl == 1
	% recover data
%     A_hat = zeros(Nr*P, Q, 'like', yp);
% 
%     for i_Q = 1 : Q
%         a_q1 = exp(-1j * 2 * pi * arg.SP * theta_hat_q(i_Q) * (0:Nr-1)).';
%         a_q2 = exp(1j * 2 * pi * eps_hat_q(i_Q) * M * (0:P-1)).';
%         A_hat(:, i_Q) = kron(a_q1, a_q2);
%     end
    a_q1    = exp(-1j * 2 * pi * arg.SP * theta_hat_q * (0:Nr-1)).';
    a_q2    = exp(1j * 2 * pi * eps_hat_q * M * (0:P-1)).';
    A_hat   = kr(a_q1, a_q2);

    % ========== Channel estimation stage 1 ==========
    if strcmpi(flags.use_chl_algo_stage_1, 'PINV')
        R_hat = A_hat \ Y_hat;
        R_hat = R_hat.';

    elseif strcmpi(flags.use_chl_algo_stage_1, 'TIKHONOV')
        R_hat = (A_hat'*A_hat + ...
                flags.chl_opt_stage_1.reg * eye(Q, 'like', yp)) \ (A_hat') ...
                * Y_hat;
        R_hat = R_hat.';

    elseif strcmpi(flags.use_chl_algo_stage_1, 'PINV_KRON')
        A_hat = kron(A_hat, eye(M, 'like', yp));
		R_hat = A_hat \ yp;
        R_hat = reshape(R_hat, M, []);

	elseif strcmpi(flags.use_chl_algo_stage_1, 'TIKHONOV_KRON')
        A_hat = kron(A_hat, eye(M, 'like', yp));
		R_hat = (A_hat'*A_hat + ...
                flags.chl_opt_stage_1.reg * eye(Q, 'like', yp)) \ (A_hat') ...
                * yp;
        R_hat = reshape(R_hat, M, []);
	else
		error("wrong algorithm selected")
    end

    % construct X_hat matrix
    X_hat = zeros(Q, M, 'like', yp);
    for i_Q = 1 : Q
        X_hat(i_Q, :) = ...
            diag(exp(-1j * 2 * pi * eps_hat_q(i_Q) * (0 : M - 1))) ...
            * R_hat(:, i_Q);
    end

    % ========== Channel estimation stage 2 ==========
	% estimate channel impulse response given training data x
    if strcmpi(flags.use_chl_algo_stage_2, 'PINV')
        h_hat_q = X_hat / (X_cir.');

	elseif strcmpi(flags.use_chl_algo_stage_2, 'MMSE')
        h_hat_q = zeros(Q, arg.Lch + 1, 'like', yp);
        for i_Q = 1 : Q
            xp_hat_q        = X_hat(i_Q, :).';
            N0_xp           = abs(xp_hat_q' * xp_hat_q - flags.chl_opt_stage_2.chlpwr) / M;
            h_hat_q(i_Q, :) = X_cir' / (X_cir * X_cir' + N0_xp * eye(M)) * xp_hat_q;
        end
		
	elseif strcmpi(flags.use_chl_algo_stage_2, 'CS')
        h_hat_q = zeros(Q, arg.Lch + 1, 'like', yp);
        chllen  = flags.chl_opt_stage_2.chllen - 1;

        % pre-compute the pseudo-inverse of small_X_cir
        n_delay = arg.Lch + 1 - chllen;
        pinv_small_X_cir = cell(n_delay, 1);
        for i_delay = 1 : n_delay
            pinv_small_X_cir{i_delay} = pinv(X_cir(:, i_delay : i_delay + chllen));
        end

        for i_Q = 1 : Q
            xp_hat_q = X_hat(i_Q, :).';

            min_fval    = inf;
            min_delay   = 0;
            min_small_h_hat = zeros(chllen + 1, 1);
            for i_delay = 1 : n_delay
                small_h_hat = pinv_small_X_cir{i_delay} * xp_hat_q;

                % evaluate 2-norm of the estimated small_h_hat
                fval = norm(xp_hat_q - X_cir(:, i_delay : i_delay + chllen) * small_h_hat)^2;
                if min_fval > fval
                    min_fval        = fval;
                    min_delay       = i_delay;
                    min_small_h_hat = small_h_hat;
                end
            end

            h_hat_q(i_Q, min_delay : min_delay + chllen) = min_small_h_hat;
        end

	else
		error("wrong algorithm selected")
    end

else
	h_hat_q = zeros(Q, arg.Lch + 1, 'like', yp);
end

end
