function y = dfo_channel_ula(x, n, eps_q, theta_q, h_q, N, Q, Nr, SP, flag_noise, add_CP_len)
%DFO_CHANNEL_ULA Generate Rx signal (in time domain)
%   x:                  input signal (in time domain)
%   n:                  additive noise
%   eps_q:              DFO
%   theta_q:            cosine value of AOA
%   h_q:                channel
%   N:                  length of x
%   Q:                  number of non-coherent source
%   Nr:                 number of RX antennas
%   SP:                 The ratio of the distance between adjacent elements
%                       in an antenna array to the wavelength of the signal
%                       (d / lambda)
%   flag_noise:         add noise or not
%   add_CP_len:         if add_CP_len <= 0, the function does not simulate 
%                       "add CP" and "remove CP" process
%

if nargin < 7
    add_CP_len = 0;
end

assert(N == length(x), "unmatch OFDM block size");
assert(Q == length(eps_q), "unmatch number of paths");
assert(Q == length(theta_q) && Q == size(h_q, 1), "unmatched Q");
assert(Nr == length(n)/N, "unmatch number of RX antennas");

% received signal
Y = zeros(N, Nr, 'like', x);

if add_CP_len > 0
    % adding CP to x, use linear convolution
    x = [x(end - add_CP_len + 1 : end) ; x];

    for m = 1 : Nr
        for i_Q = 1 : Q
            h_q_vec = h_q(i_Q, :);
            yp = conv(x.', h_q_vec);

            % add DFO
            yp = yp .* exp(1j * 2 * pi * eps_q(i_Q) * (-add_CP_len : length(yp) - add_CP_len - 1));
            
            % add AOA
            yp = yp * exp(-1j * 2 * pi * SP * theta_q(i_Q) * (m-1));

            % remove CP
            Y(:, m) = Y(:, m) + yp(add_CP_len + 1 : add_CP_len + N).';
        end
    end

else
    % not adding CP to x, use circular convolution
    for m = 1 : Nr
        for i_Q = 1 : Q
	        h_q_vec = h_q(i_Q, :);
            yp = cconv(x.', h_q_vec, N);

            % add DFO
            yp = yp .* exp(1j * 2 * pi * eps_q(i_Q) * (0 : N - 1));

            % add AOA
            yp = yp * exp(-1j * 2 * pi * SP * theta_q(i_Q) * (m-1));

	        Y(:, m) = Y(:, m) + yp.';

        end
    end

end

% vectorize Y & add noise
if flag_noise
	y = Y(:) + n;
else
	y = Y(:);
end

end
