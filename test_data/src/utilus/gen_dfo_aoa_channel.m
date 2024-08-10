function [dfoVal, aoaVal, chlVal, preset_var] = gen_dfo_aoa_channel(Q, config, preset_var)
%GEN_DFO_AOA_CHANNEL DFO & AOA channel generation
% Input arguments:
%   + Q:                # of path
%   + config:
%       - N:			# of subcarriers in one OFDM block
%       - Lcp:			CP length
%       - beta:         roll-off ratio of raised-cosine pulse
%       - fc:			carrier frequency
%       - fs:			sampling frequency
%       - mod_type:		modulation
%       - Nt:           # of Tx antenna
%       - Nr:           # of Rx antenna
%       - d:            uniform distance the antennas spaced by
%
%       - LOS:			line of sight
%       - v:			velocity of TX (km/h)
%       - DS:			delay spread
%       - pg_low:		lowest power of delay tap
%       - pg_high:		highest power of delay tap
%       - Lch:			length of channel taps
%       - doa_gap:      smallest gap between two successive angles
%       - aoa_gap:      smallest gap between two successive angles
%
%       - M:			length of subvector
%       - H:			# of center taps of RC pulse
%       - data_block:	# of data block for calculating BER

%% auxiliary variables
T           = 1. / config.fs;
% largest Q
Q_max       = max(Q);
% number of Q
n_Q         = length(Q);
% maximum Doppler shift (v: convert km/hr to m/s)
fd_q_max    = config.fc * (config.v / 3.6) / physconst('LightSpeed');
% compute minimal gap
dfo_cos_gap = config.dfo_gap;
aoa_cos_gap = config.aoa_gap;

%% main
if nargin < 3
    % 1. generate DFOs & AOAs for all possible Q's
    DFOs = gen_angle(Q_max, dfo_cos_gap, -pi, pi);
    AOAs = gen_angle(Q_max, aoa_cos_gap, 0, pi);
    
    % 2. generate normailzed gain of path (Rayleigh fading)
    CHLs = randn(Q_max, 1) + 1j * randn(Q_max, 1);
    % assume NLOS, LOS needs K (Rician/Ricean factor)
    CHLs = CHLs ./ abs(CHLs);

    % save to preset_var
    preset_var      = struct;
    preset_var.DFOs = DFOs;
    preset_var.AOAs = AOAs;
    preset_var.CHLs = CHLs;
else
    DFOs = preset_var.DFOs;
    AOAs = preset_var.AOAs;
    CHLs = preset_var.CHLs;
end

% 3. uniformly generate channel taps normalized delay
CTNDs = sort(rand(Q_max, 1));

% 4. ========== generate DFO & AOA and channel ==========
dfoVal = cell(n_Q, 1);
aoaVal = cell(n_Q, 1);
chlVal = cell(n_Q, 1);

for i_Q = 1 : n_Q
    % choose Q(i_Q) elements from Q elements
	selected_idxs = nchoosek(1 : Q_max, Q(i_Q));

    DFOq   = zeros(size(selected_idxs, 1), Q(i_Q));
    AOAq   = zeros(size(selected_idxs, 1), Q(i_Q));
    CHLq   = zeros(size(selected_idxs, 1), Q(i_Q), config.Lch + 1);

    for idx = 1 : size(selected_idxs, 1)
        selected_idx = selected_idxs(idx, :);
        
        % select DFO, AOA, CHL for this Q
        selected_dfo_q  = DFOs(selected_idx);
        selected_aoa_q  = AOAs(selected_idx);
        selected_CTND_q = CTNDs(selected_idx);
        selected_CTP_q  = get_CTP(selected_CTND_q, 1, ...
                            config.pg_low, config.pg_high);
        selected_CHL_q  = CHLs(selected_idx);
        
        % normialized channel power (path gain)
        sigma_square = 10.^(selected_CTP_q./10);
        sigma_square = sigma_square ./ sum(sigma_square);
        
        % passband channel gain (path gain)
        sigma   = sqrt(sigma_square ./ 2);
        a0_q    = sigma .* selected_CHL_q;

		% propagation delay
        tau_q = selected_CTND_q * config.DS;

		% DFO
        fd_q    = fd_q_max .* cos(selected_dfo_q);
        eps_q   = fd_q * T;

		% baseband channel gain
        a_q = a0_q .* exp(-1j * 2 * pi * (config.fc + fd_q) .* tau_q);

		% generate RC impulse of length Lch + 1 with path delay
        p_q = zeros(Q(i_Q), config.Lch + 1);
        for ii_Q = 1 : Q(i_Q)
            for i_Lch = 0 : config.Lch
                p_q(ii_Q, i_Lch + 1) = rc_pulse( ...
                    (i_Lch - (config.H - 1)/2) * T - tau_q(ii_Q), ...
                    T, ...
                    config.H, ...
                    config.beta ...
                );
            end
        end

		% channel impulse response
        h_q = a_q .* p_q;

        DFOq(idx, :)    = eps_q;
        AOAq(idx, :)    = cos(selected_aoa_q);
        CHLq(idx, :, :) = h_q;
    end

    dfoVal{i_Q} = DFOq;
    aoaVal{i_Q} = AOAq;
    chlVal{i_Q} = CHLq;
end

end

function CTP = get_CTP(CTND, CTND_max, lb, ub)
%GET_CTP Compute Channel Taps Power
%   Given the maximum value of CTND (Channel Taps Normalized Delay), the
%   CTP's are determined by the linear interpolation between lower bound
%   and upper bound.

if lb > ub;	error("lower bound > upper bound");	end
if any(CTND(:)<0) || any(CTND(:)>1);	error("CTND < 0 or CTND > 1");	end

CTND = CTND ./ CTND_max;
CTP = ub - CTND * (ub - lb);

end

function phi = gen_angle(Q, gap, a, b)
%GEN_ANGLE The function randomly generates Q angles. For any two angles,
% the absolute difference between them would be greater than gap. Note
% that it may fails to find such Q angles due to inappropriate inputs.
%

c_Q = 0;
phi = zeros(Q, 1);

c_loop = 0;
while c_Q < Q
    phi_q = a + (b - a) * rand(1);
    
    if (c_Q == 0)
        phi(1) = phi_q;
        c_Q = c_Q + 1;
    else
        too_close_flag = 0;
        % find too close angles pair
        idx = 1;
        while (idx <= c_Q)
            if (abs(cos(phi_q) - cos(phi(idx))) < gap)
                too_close_flag = 1;
                break;
            end
            idx = idx + 1;
        end

        if (too_close_flag == 0)
            phi(c_Q + 1) = phi_q;
            c_Q = c_Q + 1;
        end
    end

	% maximum loop time
	c_loop = c_loop + 1;
	if c_loop > 1000
		error("can not found!");
	end
end

end

function value = rc_pulse(t, T, H, beta)
%RC_PULSE Return the value of the Raised Cosine function at time t.
%
% t:	time
% T:	symbol time
% H:	length of center taps of RC pulse
% beta:	roll-off factor
%

h = (H-1)/2;
if abs(t) > h*T
	value = 0;
else
	if abs(t) == T/(2*beta)
    	value = (pi/4) * sinc(1/(2*beta));
	else
    	value = sinc(t/T) * cos(pi*beta*t/T) / (1 - (2*beta*t/T)^2);
	end
end

end
