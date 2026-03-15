function [hEst, delayEst, dopplerEst, navInfo] = ddChannelEstimate(rxGrid, pilotInfo)

%--------------------------------------------------------------------------
%
%   Estimates the channel in the Delay-Doppler domain using the received
%   pilot response. Also extracts navigation parameters (delay/Doppler)
%   for LEO satellite integrated communication and navigation.
%
%   Method: The transmitted pilot is a single impulse at (lp, kp).
%   After passing through an L-path channel, the received DD grid shows
%   the channel response as shifted/scaled copies of the pilot.
%   The channel taps h_i at (l_i, k_i) can be read directly from the
%   guard region of the received grid.
%
%--------------------------------------------------------------------------
% Input arguments:
%
% rxGrid            N x M received DD grid (after SFFT demodulation)
% pilotInfo         Struct from pilotPatternDD with pilot position info
%
%--------------------------------------------------------------------------
% Function returns:
%
% hEst              N x M estimated DD domain channel matrix
% delayEst          Estimated delay indices of channel paths
% dopplerEst        Estimated Doppler indices of channel paths
% navInfo           Struct with navigation-relevant measurements:
%                     .pathDelays     - delay of each detected path
%                     .pathDopplers   - Doppler shift of each path
%                     .pathGains      - complex gain of each path
%                     .snrPilot       - estimated pilot SNR
%
%--------------------------------------------------------------------------

[N, M] = size(rxGrid);
kp = pilotInfo.kp;
lp = pilotInfo.lp;
kGuard = pilotInfo.kGuard;
lGuard = pilotInfo.lGuard;

%% Extract the guard region from the received grid
% The channel response appears as copies of the pilot shifted by (l_i, k_i)
% We read the guard region centered at (lp, kp) to capture these copies.
guardRows = max(1, lp - lGuard) : min(N, lp + lGuard);
guardCols = max(1, kp - kGuard) : min(M, kp + kGuard);

guardRegion = rxGrid(guardRows, guardCols);

%% Detect channel paths via threshold
% Estimate noise power from corners of the guard region
% (where no channel response should appear if guard is large enough)
pilotEnergy = abs(rxGrid(lp, kp))^2;
noiseEst = median(abs(guardRegion(:)).^2) * 0.5;  % Robust noise estimate

% Detection threshold: paths above noise floor
thresholddB = 6;   % 6 dB above noise floor
threshold = sqrt(noiseEst * 10^(thresholddB/10));

%% Build DD channel impulse response
hEst = zeros(N, M);

% Scan the guard region for significant taps
delayEst = [];
dopplerEst = [];
pathGains = [];

for r = 1:length(guardRows)
    for c = 1:length(guardCols)
        lr = guardRows(r);
        kc = guardCols(c);
        val = rxGrid(lr, kc);

        if abs(val) > threshold
            % Channel tap detected at relative position (lr-lp, kc-kp)
            % Normalize by transmitted pilot amplitude
            hTap = val / rxGrid(lp, kp) * abs(rxGrid(lp, kp));

            % Place in DD channel matrix
            % The channel response at (l, k) means a path with
            % delay index = l - lp (mod N) and Doppler index = k - kp (mod M)
            delayShift = lr - lp;
            dopplerShift = kc - kp;

            % Map to DD channel matrix (circular shift)
            delayIdx = mod(delayShift, N) + 1;
            dopplerIdx = mod(dopplerShift, M) + 1;
            hEst(delayIdx, dopplerIdx) = val;

            delayEst = [delayEst; delayShift];
            dopplerEst = [dopplerEst; dopplerShift];
            pathGains = [pathGains; val];
        end
    end
end

%% If no paths detected, use the pilot position as single-tap channel
if isempty(delayEst)
    hEst(1, 1) = rxGrid(lp, kp);
    delayEst = 0;
    dopplerEst = 0;
    pathGains = rxGrid(lp, kp);
end

%% Fractional Doppler estimation for LoS path (parabolic interpolation)
% Find the strongest path (LoS candidate)
[~, losIdx] = max(abs(pathGains));
losRow = lp + delayEst(losIdx);   % row index in rxGrid
losCol = kp + dopplerEst(losIdx); % col index in rxGrid

% Parabolic interpolation along Doppler (column) dimension
fracDopplerShift = dopplerEst(losIdx);  % default: integer bin
if losCol > 1 && losCol < M
    alpha = abs(rxGrid(losRow, losCol - 1));
    beta  = abs(rxGrid(losRow, losCol));
    gamma = abs(rxGrid(losRow, losCol + 1));
    if (2*beta - alpha - gamma) ~= 0
        delta_k = 0.5 * (alpha - gamma) / (alpha - 2*beta + gamma);
        fracDopplerShift = dopplerEst(losIdx) + delta_k;
    end
end

% Parabolic interpolation along delay (row) dimension
fracDelayShift = delayEst(losIdx);  % default: integer bin
if losRow > 1 && losRow < N
    alpha = abs(rxGrid(losRow - 1, losCol));
    beta  = abs(rxGrid(losRow, losCol));
    gamma = abs(rxGrid(losRow + 1, losCol));
    if (2*beta - alpha - gamma) ~= 0
        delta_l = 0.5 * (alpha - gamma) / (alpha - 2*beta + gamma);
        fracDelayShift = delayEst(losIdx) + delta_l;
    end
end

%% Navigation information extraction
navInfo.pathDelays = delayEst;
navInfo.pathDopplers = dopplerEst;
navInfo.pathGains = pathGains;
navInfo.pilotPower = pilotEnergy;
navInfo.noiseEst = noiseEst;
if noiseEst > 0
    navInfo.snrPilot = 10*log10(pilotEnergy / noiseEst);
else
    navInfo.snrPilot = Inf;
end
navInfo.numPathsDetected = length(delayEst);
navInfo.losFracDoppler = fracDopplerShift;  % Fractional Doppler of LoS (bins)
navInfo.losFracDelay = fracDelayShift;      % Fractional delay of LoS (bins)

end
