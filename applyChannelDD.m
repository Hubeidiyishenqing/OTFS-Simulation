function rxDD = applyChannelDD(txDD, chInfo, scs, cpSize)

%--------------------------------------------------------------------------
%
%   Applies the multipath channel directly in the Delay-Doppler domain
%   via 2D circular convolution.
%
%   This is exact for any Doppler shift, unlike the TF-domain element-wise
%   model (channelTF .* signal) which breaks down when fd/scs is large.
%
%   The DD channel acts as: y[l,k] = sum_i h_i * x[(l-l_i) mod N, (k-k_i) mod M]
%   where l_i and k_i are the delay and Doppler bin indices of path i.
%
%--------------------------------------------------------------------------
% Input arguments:
%
% txDD              N x M transmit DD grid
% chInfo            Struct from multipathChannel with fields:
%                     .pathDelays_s    - path delays in seconds
%                     .pathDopplers_Hz - path Doppler shifts in Hz
%                     .pathGains       - path complex gains
%                     .numPaths        - number of channel paths
% scs               Subcarrier spacing (Hz)
% cpSize            Cyclic prefix ratio
%
%--------------------------------------------------------------------------
% Function returns:
%
% rxDD              N x M received DD grid after channel
%
%--------------------------------------------------------------------------

[N, M] = size(txDD);

% DD grid resolutions
Ts = (1 + cpSize) / scs;           % OFDM symbol duration with CP (s)
delta_tau = 1 / (N * scs);         % Delay resolution (s/bin)
delta_nu  = 1 / (M * Ts);          % Doppler resolution (Hz/bin)

rxDD = zeros(N, M);

for i = 1:chInfo.numPaths
    % Map physical delay/Doppler to integer bin indices
    li = round(chInfo.pathDelays_s(i) / delta_tau);     % Delay bins
    ki = round(chInfo.pathDopplers_Hz(i) / delta_nu);   % Doppler bins
    hi = chInfo.pathGains(i);

    % 2D circular convolution: shift and accumulate
    rxDD = rxDD + hi * circshift(txDD, [li, ki]);
end

end
