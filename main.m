clear;
close all;

%--------------------------------------------------------------------------
%
% This code forms a simulation of a wideband wireless communications system
% in multipath fading. It simulates: OFDM,  OTFS, coded-OFDM and coded-OTFS
% 
% The following files are required:
% dataGen.m
% multipathChannel.m
% modOFDM.m
% demodOFDM.m
% ISFFT.m
% SFFT.m
% equaliser.m
% plotGraphs.m
%
%--------------------------------------------------------------------------
%
% Author: Bradley Bates
% University of Bristol, UK
% email address: bb16177@bristol.ac.uk
% May 2020
%
% Copyright (c) 2020, Bradley Bates
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Define simulation parameters
%--------------------------------------------------------------------------
M = 16;                         % Modulation alphabet
k = log2(M);                    % Bits/symbol
cpSize = 0.07;                  % OFDM cyclic prefix size
scs = 15e3;                     % Subcarrier spacing Hz
Bw = 10e6;                      % System bandwidth Hz
ofdmSym = 14;                   % No. OFDM symbols / subframe 
EbNo = (-3:1:30)';              % Range of energy/bit to noise power ratio
velocity = 120;                 % Velocity of mobile rx relative tx km/hr
codeRate = 2/4;                 % FEC code rate used
maxIterations = 25;             % Set maximum no. of iterations for LDPC decoder
totalBits = 1e6;                % The approx. total no. of bits simulated
repeats = 1;                    % Number of simulation repetitions 


%--------------------------------------------------------------------------
%                    Initialise Simulation Components
%--------------------------------------------------------------------------

% Initialise OFDM Mod/Demod variables
numSC = pow2(ceil(log2(Bw/scs))); % Calc. nearest base-2 no. of OFDM subcarriers
cpLen = floor(cpSize * numSC);    % Calc. cyclic prefix length
numDC = (numSC - 12);             % Calc. number of data carriers

% Initialise the AWGN channel
awgnChannel = comm.AWGNChannel('NoiseMethod','Variance', 'VarianceSource','Input port');
errorRate = comm.ErrorRate('ResetInputPort',true);
errorRate1 = comm.ErrorRate('ResetInputPort',true);

% Initialise the LDPC coder/decoder
parityCheck_matrix = dvbs2ldpc(codeRate);               % Generate DVB-S.2 paraity check matrix
ldpcEncoder = comm.LDPCEncoder(parityCheck_matrix);     % Create encoder system object
ldpcDecoder = comm.LDPCDecoder(parityCheck_matrix);     % Create decoder system object
ldpcDecoder.MaximumIterationCount = maxIterations;      % Assign decoder's maximum iterations
noCodedbits = size(parityCheck_matrix,2);               % Find the Codeword length

%--------------------------------------------------------------------------
% LEO Satellite Pilot Pattern Configuration (DD domain)
%--------------------------------------------------------------------------
% For LEO: large Doppler spread requires wide Doppler guard band
pilotConfig.kp = ceil(ofdmSym/2);       % Pilot Doppler index (center)
pilotConfig.lp = ceil((Bw/scs - 12)/2); % Pilot delay index (center of data carriers)
pilotConfig.kGuard = 4;                 % Doppler guard (wide for LEO high Doppler)
pilotConfig.lGuard = 2;                 % Delay guard
pilotConfig.pilotBoostdB = 10;          % Pilot power boost (dB) - balance nav vs comm

% Create Vectors for storing error data
berOFDM = zeros(length(EbNo),3); berCOFDM = zeros(length(EbNo),3); berOTFS = zeros(length(EbNo),3); berCOTFS = zeros(length(EbNo),3);
berOTFS_pilot = zeros(length(EbNo),3); berCOTFS_pilot = zeros(length(EbNo),3);
errorStats_coded = zeros(1,3); errorStats_uncoded = zeros(1,3);

for repetition=1:repeats                                % Repeat simulation multiple times with a unqique channel for each repeat
    
    % Generate and Encode data
    [dataIn, dataBits_in, codedData_in, packetSize, numPackets, numCB] = dataGen(k,numDC,ofdmSym,totalBits,codeRate,ldpcEncoder);
    
    % Generate Rayleigh Fading Channel Impulse Response
    txSig_size = zeros((numSC+cpLen),ofdmSym);                       % Assign the size of the channel
    rayChan = multipathChannel(cpSize, scs, txSig_size, velocity);   % Create fading channel impulse response
    
    % QAM Modulator
    qamTx = qammod(dataIn,M,'InputType','bit','UnitAveragePower',true);    % Apply QAM modulation
    parallelTx = reshape(qamTx,[numDC,ofdmSym*packetSize]);                % Convert to parallel
    
    % Add nulls at index 1
    guardbandTx = [zeros(1,ofdmSym*packetSize); parallelTx];
    % Add 11 nulls around DC
    guardbandTx = [guardbandTx(1:(numDC/2),:); zeros(11,ofdmSym*packetSize); guardbandTx((numDC/2)+1:end,:)];
    
    
%--------------------------------------------------------------------------
%                       OFDM BER Calculation
%--------------------------------------------------------------------------
    
    % Calculate SNR
    snr = EbNo + 10*log10(codeRate*k) + 10*log10(numDC/((numSC)));
    
    % Multicarrier Modulation
    frameBuffer = guardbandTx;          % Create a 'buffer' so subframes can be individually modulated
    txframeBuffer = [];                 % Initilaise matrix
    for w = 1:packetSize
        ofdmTx = modOFDM(frameBuffer(:,1:ofdmSym),numSC,cpLen,ofdmSym);    % Apply OFDM modulation to a subframe of data
        frameBuffer(:, 1:ofdmSym) = [];                                    % Delete modulated data from frameBuffer
        txframeBuffer = [txframeBuffer;ofdmTx];                            % Add modulated subframe to transmit buffer
    end
    
    
    % Loop through different values of EbNo
    for m = 1:length(EbNo)
        % Loop through the of packets to be transmitted
        for j = 1:numPackets
            rxframeBuffer = [];                 % Initialise matrix
            
            % Transmit each subframe individually
            for u = 1:packetSize
                
                % Remove next subframe from the transmit buffer
                txSig = txframeBuffer( ((u-1)*numel(ofdmTx)+1) : u*numel(ofdmTx) );
                
                % Apply Channel to input signal
                fadedSig = zeros(size(txSig));                    % Pre-allocate vector size
                for i = 1:size(txSig,1)                           % Perform elementwise...
                    for j = 1:size(txSig,2)                       % ...matrix multiplication
                        fadedSig(i,j) = txSig(i,j).*rayChan(i,j);
                    end
                end
                
                % AWGN Channel
                release(awgnChannel);
                powerDB = 10*log10(var(fadedSig));            % Calculate Tx signal power
                noiseVar = 10.^(0.1*(powerDB-snr(m)));        % Calculate the noise variance
                rxSig = awgnChannel(fadedSig,noiseVar);       % Pass the signal through a noisy channel
                
                % Equalisation
                eqSig = equaliser(rxSig,fadedSig,txSig,ofdmSym);
                
                % Demodulation
                rxSubframe = demodOFDM(eqSig,cpLen,ofdmSym);     % Apply OFDM demodulation
                rxframeBuffer = [rxframeBuffer';rxSubframe']';         % Store demodulated subframe in rx buffer
            end
            % Remove all null carriers
            parallelRx = rxframeBuffer;
            parallelRx((numDC/2)+1:(numDC/2)+11, :) = [];     % Remove nulls around the DC input
            parallelRx(1:1, :) = [];                          % Remove nulls at index 1
            qamRx = reshape(parallelRx,[numel(parallelRx),1]);% Convert to serial
            
            % Uncoded demodulation of entire packet
            dataOut = qamdemod(qamRx,M,'OutputType','bit','UnitAveragePower',true);% Apply QAM demodulation
            codedData_out = randdeintrlv(dataOut,4831);                            % De-interleave data
            codedData_out(numel(codedData_in)+1:end) = [];                         % Remove pad bits
            errorStats_uncoded = errorRate(codedData_in,codedData_out,0);          % Collect error statistics
          
            % Coded demodulation of entire packet
            powerDB = 10*log10(var(qamRx));                                   % Calculate Rx signal power
            noiseVar = 10.^(0.1*(powerDB-(EbNo(m) + 10*log10(codeRate*k) - 10*log10(sqrt(numDC)))));            % Calculate the noise variance
            dataOut = qamdemod(qamRx,M,'OutputType', 'approxllr','UnitAveragePower',true,'NoiseVariance',noiseVar);% Apply QAM demodulation
            codedData_out1 = randdeintrlv(dataOut,4831);                      % De-interleave data
            codedData_out1(numel(codedData_in)+1:end) = [];                   % Remove pad bits
            
            % Decode individual code blocks
            dataBits_out = [];                                                % Initialise matrix
            dataOut_buffer = codedData_out1;
            for q = 1:numCB
                dataBits_out = [dataBits_out;ldpcDecoder(dataOut_buffer(1:noCodedbits))]; % Decode data & add it to the data bits out matrix
                dataOut_buffer(1:noCodedbits) = [];                                       % Delete decoded data from buffer
            end
            dataBits_out = double(dataBits_out);                              % Convert to a double compatible w/ errorStats
            errorStats_coded = errorRate1(dataBits_in,dataBits_out,0);        % Collect error statistics
            
        end
        berOFDM(m,:) = errorStats_uncoded;                                  % Save uncoded BER data
        berCOFDM(m,:) = errorStats_coded;                                   % Save  coded BER data
        errorStats_uncoded = errorRate(codedData_in,codedData_out,1);       % Reset the error rate calculator
        errorStats_coded = errorRate1(dataBits_in,dataBits_out,1);          % Reset the error rate calculator
        
    end
    
    
%--------------------------------------------------------------------------
%                       OTFS BER Calculation
%--------------------------------------------------------------------------
    
    % Calculate SNR
    snr = EbNo + 10*log10(codeRate*k) + 10*log10(numDC/((numSC))) + 10*log10(sqrt(ofdmSym));
    
    % Multicarrier Modulation
    frameBuffer = guardbandTx;          % Create a 'buffer' so subframes can be individually modulated
    txframeBuffer = [];                 % Initilaise matrix
    for w = 1:packetSize
        otfsTx = ISFFT(frameBuffer(:,1:ofdmSym));       % Apply OTFS modulation to a subframe of data
        ofdmTx = modOFDM(otfsTx,numSC,cpLen,ofdmSym);    % Apply OFDM modulation
        frameBuffer(:, 1:ofdmSym) = [];                  % Delete modulated data from frameBuffer
        txframeBuffer = [txframeBuffer;ofdmTx];          % Add modulated subframe to transmit buffer
    end
    
    % Loop through different values of EbNo
    for m = 1:length(EbNo)
        % Loop through the of packets to be transmitted
        for j = 1:numPackets
            rxframeBuffer = [];                 % Initialise matrix
            
            % Transmit each subframe individually
            for u = 1:packetSize
                
                % Remove next subframe from the transmit buffer
                txSig = txframeBuffer( ((u-1)*numel(ofdmTx)+1) : u*numel(ofdmTx) );
                
                % Apply Channel to input signal
                fadedSig = zeros(size(txSig));                    % Pre-allocate vector size
                for i = 1:size(txSig,1)                           % Perform elementwise...
                    for j = 1:size(txSig,2)                       % ...matrix multiplication
                        fadedSig(i,j) = txSig(i,j).*rayChan(i,j);
                    end
                end
                
                % AWGN Channel
                release(awgnChannel);
                powerDB = 10*log10(var(fadedSig));            % Calculate Tx signal power
                noiseVar = 10.^(0.1*(powerDB-snr(m)));        % Calculate the noise variance
                rxSig = awgnChannel(fadedSig,noiseVar);       % Pass the signal through a noisy channel
                
                % Equalisation
                eqSig = equaliser(rxSig,fadedSig,txSig,ofdmSym);
                
                % Demodulation
                otfsRx = demodOFDM(eqSig,cpLen,ofdmSym);     % Apply OFDM demodulation
                rxSubframe = SFFT(otfsRx);                      % Apply OTFS demodulation
                rxframeBuffer = [rxframeBuffer';rxSubframe']';     % Store demodulated subframe in rx buffer
            end
            % Remove all null carriers
            parallelRx = rxframeBuffer;
            parallelRx((numDC/2)+1:(numDC/2)+11, :) = [];         % Remove nulls around the DC input
            parallelRx(1:1, :) = [];                              % Remove nulls at index 1
            qamRx = reshape(parallelRx,[numel(parallelRx),1]);    % Convert to serial
            
            % Uncoded demodulation of entire packet
            dataOut = qamdemod(qamRx,M,'OutputType','bit','UnitAveragePower',true);% Apply QAM demodulation
            codedData_out = randdeintrlv(dataOut,4831);                            % De-interleave data
            codedData_out(numel(codedData_in)+1:end) = [];                         % Remove pad bits
            errorStats_uncoded = errorRate(codedData_in,codedData_out,0);          % Collect error statistics
            

            % Coded demodulation of entire packet
            powerDB = 10*log10(var(qamRx));                                   % Calculate Rx signal power
            noiseVar = 10.^(0.1*(powerDB-(EbNo(m) + 10*log10(codeRate*k) - 10*log10(sqrt(numDC)))));            % Calculate the noise variance
            dataOut = qamdemod(qamRx,M,'OutputType', 'approxllr','UnitAveragePower',true,'NoiseVariance',noiseVar);% Apply QAM demodulation
            codedData_out1 = randdeintrlv(dataOut,4831);                      % De-interleave data
            codedData_out1(numel(codedData_in)+1:end) = [];                   % Remove pad bits
            
            % Decode individual code blocks
            dataBits_out = [];                                                % Initialise matrix
            dataOut_buffer = codedData_out1;
            for q = 1:numCB
                dataBits_out = [dataBits_out;ldpcDecoder(dataOut_buffer(1:noCodedbits))]; % Decode data & add it to the data bits out matrix
                dataOut_buffer(1:noCodedbits) = [];                                       % Delete decoded data from buffer
            end
            dataBits_out = double(dataBits_out);                              % Convert to a double compatible w/ errorStats
            errorStats_coded = errorRate1(dataBits_in,dataBits_out,0);     % Collect error statistics
            
        end
        berOTFS(m,:) = errorStats_uncoded;                                  % Save uncoded BER data
        berCOTFS(m,:) = errorStats_coded;                                   % Save coded BER data
        errorStats_uncoded = errorRate(codedData_in,codedData_out,1);       % Reset the error rate calculator
        errorStats_coded = errorRate1(dataBits_in,dataBits_out,1);          % Reset the error rate calculator

    end


%--------------------------------------------------------------------------
%              OTFS with DD-Domain Pilot Pattern (LEO Satellite)
%--------------------------------------------------------------------------

    % Calculate SNR (same as OTFS)
    snr = EbNo + 10*log10(codeRate*k) + 10*log10(numDC/((numSC))) + 10*log10(sqrt(ofdmSym));

    % --- Transmitter: Place pilot + data on DD grid ---
    % For each subframe, create DD grid with pilot pattern
    frameBuffer_coded = guardbandTx;    % Use same guard-banded QAM data
    txframeBuffer_pilot = [];

    % Store pilot info for all subframes (same pattern each subframe)
    pilotConfig.lp = ceil(size(guardbandTx,1)/2);  % Update to match actual grid rows
    pilotConfig.kp = ceil(ofdmSym/2);

    for w = 1:packetSize
        subframe = frameBuffer_coded(:, 1:ofdmSym);
        frameBuffer_coded(:, 1:ofdmSym) = [];

        % Extract data symbols from subframe (column-major order)
        dataSyms = subframe(:);

        % Create DD grid with pilot pattern
        [ddGrid, dataIdx, pilotIdx, guardIdx, pilotInfoTx] = ...
            pilotPatternDD(dataSyms, size(subframe,1), ofdmSym, pilotConfig);

        % OTFS modulation: ISFFT (DD -> TF) then OFDM modulation
        otfsTx = ISFFT(ddGrid);
        ofdmTx = modOFDM(otfsTx, numSC, cpLen, ofdmSym);
        txframeBuffer_pilot = [txframeBuffer_pilot; ofdmTx];
    end

    % --- Channel + Receiver ---
    for m = 1:length(EbNo)
        for j = 1:numPackets
            rxframeBuffer = [];

            for u = 1:packetSize
                % Extract subframe
                txSig = txframeBuffer_pilot( ((u-1)*numel(ofdmTx)+1) : u*numel(ofdmTx) );

                % Apply multipath fading channel
                fadedSig = zeros(size(txSig));
                for i = 1:size(txSig,1)
                    for jj = 1:size(txSig,2)
                        fadedSig(i,jj) = txSig(i,jj).*rayChan(i,jj);
                    end
                end

                % AWGN Channel
                release(awgnChannel);
                powerDB = 10*log10(var(fadedSig));
                noiseVar = 10.^(0.1*(powerDB-snr(m)));
                rxSig = awgnChannel(fadedSig, noiseVar);

                % OFDM demodulation -> SFFT (TF -> DD)
                % Use time-domain equalizer first, then SFFT
                eqSig = equaliser(rxSig, fadedSig, txSig, ofdmSym);
                otfsRx = demodOFDM(eqSig, cpLen, ofdmSym);
                rxDD = SFFT(otfsRx);

                % DD-domain channel estimation using pilot
                [hEst, delayEst, dopplerEst, navInfo] = ddChannelEstimate(rxDD, pilotInfoTx);

                % DD-domain equalization
                [eqDD, eqDataSyms] = otfsEqualiserDD(rxDD, hEst, dataIdx, pilotInfoTx);

                % Reconstruct subframe from equalized data at data positions
                rxSubframe = zeros(size(subframe));
                if length(eqDataSyms) >= length(dataIdx)
                    rxSubframe(dataIdx) = eqDataSyms(1:length(dataIdx));
                end

                rxframeBuffer = [rxframeBuffer'; rxSubframe']';
            end

            % Remove guard band nulls (same as original OTFS path)
            parallelRx = rxframeBuffer;
            parallelRx((numDC/2)+1:(numDC/2)+11, :) = [];
            parallelRx(1:1, :) = [];
            qamRx = reshape(parallelRx, [numel(parallelRx), 1]);

            % Uncoded demodulation
            dataOut = qamdemod(qamRx, M, 'OutputType', 'bit', 'UnitAveragePower', true);
            codedData_out = randdeintrlv(dataOut, 4831);
            codedData_out(numel(codedData_in)+1:end) = [];
            errorStats_uncoded = errorRate(codedData_in, codedData_out, 0);

            % Coded demodulation
            powerDB = 10*log10(var(qamRx));
            noiseVar = 10.^(0.1*(powerDB-(EbNo(m) + 10*log10(codeRate*k) - 10*log10(sqrt(numDC)))));
            dataOut = qamdemod(qamRx, M, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', noiseVar);
            codedData_out1 = randdeintrlv(dataOut, 4831);
            codedData_out1(numel(codedData_in)+1:end) = [];

            dataBits_out = [];
            dataOut_buffer = codedData_out1;
            for q = 1:numCB
                dataBits_out = [dataBits_out; ldpcDecoder(dataOut_buffer(1:noCodedbits))];
                dataOut_buffer(1:noCodedbits) = [];
            end
            dataBits_out = double(dataBits_out);
            errorStats_coded = errorRate1(dataBits_in, dataBits_out, 0);
        end

        berOTFS_pilot(m,:) = errorStats_uncoded;
        berCOTFS_pilot(m,:) = errorStats_coded;
        errorStats_uncoded = errorRate(codedData_in, codedData_out, 1);
        errorStats_coded = errorRate1(dataBits_in, dataBits_out, 1);
    end

    % Print navigation info from last subframe
    fprintf('\n--- LEO Navigation Info (from DD-domain pilot) ---\n');
    fprintf('Detected %d channel paths\n', navInfo.numPathsDetected);
    fprintf('Pilot SNR: %.1f dB\n', navInfo.snrPilot);
    fprintf('Pilot overhead: %.1f%%\n', pilotInfoTx.overheadPercent);
    fprintf('Effective pilot boost: %.1f dB\n', pilotInfoTx.effectivePilotBoostdB);
    for pp = 1:length(navInfo.pathDelays)
        fprintf('  Path %d: delay=%+d bins, Doppler=%+d bins, |gain|=%.3f\n', ...
            pp, navInfo.pathDelays(pp), navInfo.pathDopplers(pp), abs(navInfo.pathGains(pp)));
    end
    fprintf('---------------------------------------------------\n');

end

%--------------------------------------------------------------------------
%                           Figures
%-------------------------------------------------------------------------- 

% Plot BER / EbNo curves (including OTFS-Pilot results)
plotGraphs(berOFDM, berCOFDM, berOTFS, berCOTFS, M, numSC, EbNo, berOTFS_pilot, berCOTFS_pilot);

% Plot DD-domain pilot pattern visualization
figure;
ddVis = zeros(size(guardbandTx,1), ofdmSym);
ddVis(dataIdx) = 0.3;                  % Data positions (gray)
ddVis(guardIdx) = 0;                   % Guard band (black)
ddVis(pilotIdx) = 1;                   % Pilot (bright)
imagesc(ddVis);
colormap(hot);
colorbar;
title('Delay-Doppler Grid: Pilot Pattern with Guard Bands');
xlabel('Doppler Index (OFDM symbols)');
ylabel('Delay Index (Subcarriers)');
hold on;
% Draw guard band rectangle
rectangle('Position', [pilotConfig.kp-pilotConfig.kGuard-0.5, ...
    pilotConfig.lp-pilotConfig.lGuard-0.5, ...
    2*pilotConfig.kGuard+1, 2*pilotConfig.lGuard+1], ...
    'EdgeColor', 'c', 'LineWidth', 2, 'LineStyle', '--');
legend('Guard Band Boundary');
hold off;

